#!/usr/bin/env bash
# run_differential.sh -- the differential pipeline, end to end.
#
#   GTEx tissue (recount3)  ->  tissue reference profile
#   cohorts                 ->  deviation from it
#   spectra / condition / consensus / library on the deviation, bp axis
#   crest genes per characteristic peak; gene-level LOCO next to the spectral LOCO
#
# THREE LIBRARIES, NOT ONE
# ------------------------
# Subtracting a recount3/GTEx profile from a matrix quantified by another
# pipeline subtracts pipeline as much as biology. So:
#   primary      recount3-only cohorts (GTEx + SRPs recount3 holds). The claim.
#   sensitivity  GTEx + the GEO-only cohorts. Shows what the pipeline
#                difference does; never the headline number.
#   combined     everything, built only with TSF_RUN_COMBINED=1 and only worth
#                reading if the LOCO of primary and sensitivity say technical
#                origin does not dominate the classification.
# Each library has its own interim tree (its own shared mask) and results tree.
#
# RESUMABLE BY CONTENT, NOT BY EXISTENCE
# --------------------------------------
# Every reusable artefact carries a sidecar <artefact>.inputs with a digest of
# everything that went into it: dataset list, GTEx project, tissue, GSE list,
# vocabulary, axis, bin size, profile digest, and the md5 of the EFFECTIVE
# project configuration (config/project.R as R loads it, every TSF_* override
# applied -- scripts/config_digest.R). The artefact is rebuilt when that digest
# differs; tsf stages get --force in that case, since they reuse their own
# files otherwise; directories are rebuilt into .tmp and swapped in with the
# previous version kept as .bak.<timestamp>. Re-run the same command after a
# failure and it continues.
#
# Flags (every one has a default; the environment is only a fallback):
#   --root DIR          repo checkout                 (default: directory of this script/..)
#   --geo-dir DIR       raw inputs                    (required)
#   --interim-dir DIR   interim trees, NEW location   (required)
#   --results-dir DIR   results trees, NEW location   (required)
#   --gse LIST          GEO cohorts                   (GSE135251,GSE130970,GSE162694,GSE276114,GSE142530)
#   --gtex PROJECT      recount3 GTEx project         (LIVER)
#   --tissue LABEL      tissue label                  (liver)
#   --vocab ID          vocabulary; must hold Control_external_study (liver_fibrosis)
#   --bin-size N        bp bin width                  (100000)
#   --workers N         cores for the library builder (4)
#   --condition-b N     permutations for the condition test (config default)
#   --maxt-b N          permutations for per-sample maxT     (config default)
#   --combined          also build the combined library
#   --only LIB          primary | sensitivity_geo | combined: build just that one
#   --skip-tests        do not run `make test` in step 0 (you ran it by hand; logged)
#   -h, --help
#
# Usage:
#   bash scripts/run_differential.sh --geo-dir /d/geo --interim-dir /d/interim_diff \
#        --results-dir /d/results_diff --workers 8 2>&1 | tee /d/results_diff/run.log
#
set -euo pipefail

usage() { sed -n '2,/^set -euo/p' "$0" | grep '^#' | sed 's/^# \{0,1\}//'; }
TSF_ROOT="${TSF_ROOT:-$(cd "$(dirname "$0")/.." && pwd)}"
TSF_GSE="${TSF_GSE:-GSE135251,GSE130970,GSE162694,GSE276114,GSE142530}"
TSF_GTEX_TISSUE="${TSF_GTEX_TISSUE:-LIVER}"
TSF_TISSUE="${TSF_TISSUE:-liver}"
TSF_VOCAB="${TSF_VOCAB:-liver_fibrosis}"
TSF_BIN_SIZE="${TSF_BIN_SIZE:-100000}"
N_WORKERS="${N_WORKERS:-4}"
TSF_RUN_COMBINED="${TSF_RUN_COMBINED:-0}"
ONLY=""; SKIP_TESTS=0
while [ $# -gt 0 ]; do
  case "$1" in
    --root)        TSF_ROOT="$2"; shift 2 ;;
    --geo-dir)     TSF_GEO_DIR="$2"; shift 2 ;;
    --interim-dir) TSF_INTERIM_DIR="$2"; shift 2 ;;
    --results-dir) TSF_RESULTS_DIR="$2"; shift 2 ;;
    --gse)         TSF_GSE="$2"; shift 2 ;;
    --gtex)        TSF_GTEX_TISSUE="$2"; shift 2 ;;
    --tissue)      TSF_TISSUE="$2"; shift 2 ;;
    --vocab)       TSF_VOCAB="$2"; shift 2 ;;
    --bin-size)    TSF_BIN_SIZE="$2"; shift 2 ;;
    --workers)     N_WORKERS="$2"; shift 2 ;;
    --condition-b) export TSF_CONDITION_B="$2"; shift 2 ;;
    --maxt-b)      export TSF_MAXT_B="$2"; shift 2 ;;
    --combined)    TSF_RUN_COMBINED=1; shift ;;
    --only)        ONLY="$2"; shift 2 ;;
    --skip-tests)  SKIP_TESTS=1; shift ;;
    -h|--help)     usage; exit 0 ;;
    *) echo "unknown flag: $1" >&2; usage >&2; exit 1 ;;
  esac
done
# relative paths become absolute: tsf and the R scripts record them in
# markers and manifests, and a relative one breaks the md5 checks later
abs() { case "$1" in /*) printf '%s' "$1" ;; *) printf '%s/%s' "$PWD" "${1#./}" ;; esac; }
[ -z "${TSF_GEO_DIR:-}" ]     || TSF_GEO_DIR=$(abs "$TSF_GEO_DIR")
[ -z "${TSF_INTERIM_DIR:-}" ] || TSF_INTERIM_DIR=$(abs "$TSF_INTERIM_DIR")
[ -z "${TSF_RESULTS_DIR:-}" ] || TSF_RESULTS_DIR=$(abs "$TSF_RESULTS_DIR")
: "${TSF_GEO_DIR:?--geo-dir is required}"
: "${TSF_INTERIM_DIR:?--interim-dir is required}"
: "${TSF_RESULTS_DIR:?--results-dir is required}"
# tsf and the R scripts read these from the environment; the flags are the
# only interface the user needs to touch
export TSF_ROOT TSF_GEO_DIR TSF_INTERIM_DIR TSF_RESULTS_DIR
cd "$TSF_ROOT"
ROOT_INTERIM="$TSF_INTERIM_DIR"; ROOT_RESULTS="$TSF_RESULTS_DIR"
mkdir -p "$TSF_GEO_DIR" "$ROOT_INTERIM" "$ROOT_RESULTS"
log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
step() { log "=== $* ==="; }
md5()  { if command -v md5sum >/dev/null 2>&1; then md5sum "$1" | awk '{print $1}'; else md5 -q "$1"; fi; }

# ---------------------------------------------------------------- inputs digest
# The EFFECTIVE configuration: config/project.R loaded by R with every TSF_*
# environment override applied (estimator, primary scheme, multitaper NW/K,
# bin aggregate and coverage, annotation format, stability criterion, ...),
# minus the path fields. Grepping the file could not see the overrides.
CONFIG_MD5=$(Rscript scripts/config_digest.R 2>/dev/null | tail -1 | tr -d ' ')
[ -n "$CONFIG_MD5" ] || { echo "config_digest.R produced nothing"; exit 1; }
log "effective config digest $CONFIG_MD5"
base_digest() {  # $1 = dataset list
  printf 'datasets=%s|gtex=%s|tissue=%s|gse=%s|vocab=%s|axis=bp|bin=%s|config=%s' \
    "$1" "$TSF_GTEX_TISSUE" "$TSF_TISSUE" "$TSF_GSE" "$TSF_VOCAB" "$TSF_BIN_SIZE" "$CONFIG_MD5"
}
# fresh ARTEFACT DIGEST : 0 if artefact exists and its .inputs equals DIGEST
fresh() { [ -e "$1" ] && [ -f "$1.inputs" ] && [ "$(cat "$1.inputs")" = "$2" ]; }
stamp() { printf '%s' "$2" > "$1.inputs"; }
# ensure ARTEFACT DIGEST CMD... : run CMD unless fresh, then stamp
ensure() {
  local art="$1" dig="$2"; shift 2
  if fresh "$art" "$dig"; then log "reuse $(basename "$art") (inputs unchanged)"; return; fi
  [ -e "$art" ] && log "rebuild $(basename "$art"): inputs changed"
  "$@"; stamp "$art" "$dig"
}
# stale ARTEFACT DIGEST : 0 if the artefact exists but its inputs changed.
# Used to add --force to tsf stages that otherwise reuse their own files
# ("reusing existing maxT"): the marker would then say "new B" while some
# products still came from the old one. --force only when something changed.
stale() { [ -e "$1" ] && ! fresh "$1" "$2"; }
# ensure_dir_atomic DIR DIGEST OUTFLAG CMD... : build into DIR.tmp via
# "CMD... OUTFLAG DIR.tmp", then swap in; the previous DIR is kept as
# DIR.bak.<timestamp>. Old signatures of a condition the new run no longer
# produces cannot survive in the new directory.
ensure_dir_atomic() {
  local dir="$1" dig="$2" outflag="$3"; shift 3
  if fresh "$dir" "$dig"; then log "reuse $(basename "$dir") (inputs unchanged)"; return; fi
  rm -rf "$dir.tmp"; mkdir -p "$dir.tmp"
  "$@" "$outflag" "$dir.tmp"
  if [ -d "$dir" ]; then
    local bak="$dir.bak.$(date +%Y%m%d_%H%M%S)"; mv "$dir" "$bak"; log "previous $(basename "$dir") kept as $(basename "$bak")"
  fi
  mv "$dir.tmp" "$dir"; stamp "$dir" "$dig"
}

# ---------------------------------------------------------------- 0. sanity
step "0 tests and recount3 config validation"
TESTS_DIG="tests=$(git rev-parse HEAD 2>/dev/null || echo nogit)"
if [ "$SKIP_TESTS" = "1" ]; then
  log "tests SKIPPED by --skip-tests (HEAD $(git rev-parse --short HEAD 2>/dev/null || echo nogit)); run make test by hand"
  printf 'skipped by flag at %s, HEAD %s\n' "$(date)" "$(git rev-parse HEAD 2>/dev/null || echo nogit)" >> "$ROOT_RESULTS/tests_skipped.log"
else
  # the tests must see a CLEAN environment: test_labels.R checks what happens
  # when no TSF_* path is set, and this script has just exported them. Output
  # goes to tests.log, and its tail is printed on failure.
  ensure "$ROOT_RESULTS/.tests_ok" "$TESTS_DIG" bash -c \
    'cd "$2" && env -u TSF_GEO_DIR -u TSF_INTERIM_DIR -u TSF_RESULTS_DIR -u TSF_CONFIG make test > "$1" 2>&1 && touch "$0" || { echo "tests failed; last lines of $1:"; tail -8 "$1"; exit 1; }' \
    "$ROOT_RESULTS/.tests_ok" "$ROOT_RESULTS/tests.log" "$TSF_ROOT"
fi
# the recount3 configs THIS run uses are validated once the dataset list is
# known (step 2); other tissues' R3_*.R in the repo are not this run's concern

# ---------------------------------------------------------- 1. which sources
step "1 resolve GEO -> SRP -> recount3"
SOURCES="config/recount3_sources.tsv"
ensure "$SOURCES" "gse=$TSF_GSE" Rscript scripts/resolve_recount3.R --gse "$TSF_GSE" --out "$SOURCES"
SRPS=$(awk -F'\t' 'NR>1 && $3=="TRUE" {print $2}' "$SOURCES" | paste -sd, -)
GEO_ONLY=$(awk -F'\t' 'NR>1 && $3!="TRUE" {print $1}' "$SOURCES" | paste -sd, -)
log "recount3 cohorts: ${SRPS:-none}   GEO-only cohorts: ${GEO_ONLY:-none}"

# ------------------------------------------------------------- 2. recount3
GTEX_ID="R3_${TSF_GTEX_TISSUE}"
step "2 fetch recount3 ($GTEX_ID${SRPS:+, $SRPS})"
ensure "$TSF_GEO_DIR/${GTEX_ID}_reads.tsv.gz" "gtex=$TSF_GTEX_TISSUE|tissue=$TSF_TISSUE|vocab=$TSF_VOCAB|fetch=v2" \
  Rscript scripts/recount3_fetch.R --projects "$TSF_GTEX_TISSUE" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
R3_COHORTS=""
if [ -n "$SRPS" ]; then
  for s in ${SRPS//,/ }; do
    ensure "$TSF_GEO_DIR/R3_${s}_reads.tsv.gz" "srp=$s|tissue=$TSF_TISSUE|vocab=$TSF_VOCAB|fetch=v2" \
      Rscript scripts/recount3_fetch.R --projects "$s" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
    if grep -q '^\s*# list(id = "biopsy_fibrosis_stage"' "config/datasets/R3_$s.R"; then
      # recount3 carried no usable sample attributes: take the labels from the
      # GEO series matrix of the same study (joined on SRX) and copy the GEO
      # config's rules verbatim. Only if that fails is a hand edit needed.
      gse=$(awk -F'\t' -v srp="$s" 'NR>1 && $2==srp {print $1}' "$SOURCES" | head -1)
      if [ -n "$gse" ] && [ -f "config/datasets/$gse.R" ]; then
        ./tsf fetch "$gse" >/dev/null 2>&1 || true
        if Rscript scripts/recount3_join_geo.R --dataset "R3_$s" --geo "$gse" --geo-dir "$TSF_GEO_DIR"; then
          log "R3_$s labelled from $gse (series matrix joined on SRX; rules copied from config/datasets/$gse.R)"
        else
          log "STOP: could not join R3_$s to $gse. Edit config/datasets/R3_$s.R by hand, then re-run."; exit 2
        fi
      else
        log "STOP: config/datasets/R3_$s.R needs condition_rules and no GEO config is known for $s. Edit it, then re-run."; exit 2
      fi
    fi
    R3_COHORTS="${R3_COHORTS:+$R3_COHORTS,}R3_$s"
  done
fi
step "2b validate the recount3 configs this run uses"
Rscript scripts/validate_recount3_configs.R --datasets "$GTEX_ID${R3_COHORTS:+,$R3_COHORTS}" --vocabulary "$TSF_VOCAB" || {
  log "STOP: review the configs above, then:"
  log "      Rscript scripts/validate_recount3_configs.R --datasets $GTEX_ID${R3_COHORTS:+,$R3_COHORTS} --migrate --vocabulary $TSF_VOCAB"
  log "      (header fields only; condition_rules are never touched; .bak kept)"
  exit 2
}

# ---------------------------------------------------------- 3. GEO inputs
if [ -n "$GEO_ONLY" ]; then
  step "3 GEO fetch/check for $GEO_ONLY"
  ./tsf fetch ${GEO_ONLY//,/ } || true
  ./tsf check ${GEO_ONLY//,/ }
fi

# --------------------------------------------------------------- libraries
# run_library NAME DATASETS  -- one interim tree, one mask, one profile, one results tree
run_library() {
  local NAME="$1" COHORTS="$2"
  local DS_LIST=${COHORTS//,/ }
  export TSF_INTERIM_DIR="$ROOT_INTERIM/$NAME" TSF_RESULTS_DIR="$ROOT_RESULTS/$NAME"
  mkdir -p "$TSF_INTERIM_DIR" "$TSF_RESULTS_DIR"
  local DIG; DIG=$(base_digest "$COHORTS")
  local MASK="$TSF_INTERIM_DIR/shared_gene_mask.tsv"
  local PROFILE="$TSF_INTERIM_DIR/reference_profile_${TSF_TISSUE}.tsv"
  step "[$NAME] datasets: $COHORTS"
  printf '%s\n' "$DIG" > "$TSF_RESULTS_DIR/inputs_digest.txt"

  step "[$NAME] ingest pass 1 (per-cohort expression filter)"
  for d in $DS_LIST; do
    ensure "$TSF_INTERIM_DIR/$d/retained_genes.tsv" "$DIG|pass=1|ds=$d" \
      ./tsf ingest "$d" --grid-axis bp --bin-size "$TSF_BIN_SIZE"
  done
  step "[$NAME] shared gene mask over exactly these cohorts"
  ensure "$MASK" "$DIG|pass=mask" \
    Rscript scripts/shared_gene_mask.R --interim-dir "$TSF_INTERIM_DIR" --datasets "$COHORTS" --out "$MASK"
  local MASK_MD5; MASK_MD5=$(md5 "$MASK")
  step "[$NAME] ingest pass 2 with the shared mask"
  for d in $DS_LIST; do
    ensure "$TSF_INTERIM_DIR/$d/.masked" "$DIG|pass=2|mask=$MASK_MD5|ds=$d" \
      bash -c './tsf ingest "$1" --grid-axis bp --bin-size "$2" --gene-mask "$3" --force && touch "$4"' _ \
        "$d" "$TSF_BIN_SIZE" "$MASK" "$TSF_INTERIM_DIR/$d/.masked"
  done

  step "[$NAME] tissue reference profile from $GTEX_ID"
  ensure "$PROFILE" "$DIG|mask=$MASK_MD5|ref=$GTEX_ID" \
    Rscript scripts/build_tissue_reference.R --datasets "$GTEX_ID" --tissue "$TSF_TISSUE" --out "$PROFILE"
  local PROF_MD5; PROF_MD5=$(md5 "$PROFILE")
  step "[$NAME] every cohort -> deviation from the profile (md5 ${PROF_MD5:0:8}); GTEx included, its deviations are the healthy class"
  # content-checked inside (md5 of expression.tsv vs marker): safe to call every time
  Rscript scripts/apply_reference_profile.R --datasets "$COHORTS" --profile "$PROFILE"

  step "[$NAME] spectra ... compare on the deviation (null = all, bp axis)"
  local STAGES_DIG="$DIG|mask=$MASK_MD5|profile=$PROF_MD5|B=${TSF_CONDITION_B:-cfg}/${TSF_MAXT_B:-cfg}"
  local FORCE=""
  if stale "$TSF_RESULTS_DIR/.stages_done" "$STAGES_DIG"; then
    FORCE="--force"; log "[$NAME] methodological inputs changed since the last run: stages will recompute (--force)"
  fi
  ensure "$TSF_RESULTS_DIR/.stages_done" "$STAGES_DIG" \
    bash -c './tsf run $1 --from spectra --to compare --grid-axis bp --bin-size "$2" --gene-mask "$3" --stage-order F0,F1,F2,F3,F4 $5 && touch "$4"' _ \
      "$DS_LIST" "$TSF_BIN_SIZE" "$MASK" "$TSF_RESULTS_DIR/.stages_done" "$FORCE"
  step "[$NAME] fingerprint library + out-of-cohort validation (profile stored in reference.rds)"
  local REF_DIG="$DIG|mask=$MASK_MD5|profile=$PROF_MD5"
  if stale "$TSF_RESULTS_DIR/reference/reference.rds" "$REF_DIG"; then
    local rbak="$TSF_RESULTS_DIR/reference.bak.$(date +%Y%m%d_%H%M%S)"; mv "$TSF_RESULTS_DIR/reference" "$rbak"
    log "[$NAME] previous reference kept as $(basename "$rbak")"
  fi
  ensure "$TSF_RESULTS_DIR/reference/reference.rds" "$REF_DIG" \
    ./tsf reference $DS_LIST --grid-axis bp --bin-size "$TSF_BIN_SIZE" --gene-mask "$MASK"

  step "[$NAME] condition library (cross-cohort meta-analysis), built atomically"
  ensure_dir_atomic "$TSF_RESULTS_DIR/condition_library" "$REF_DIG" --out-dir \
    Rscript scripts/build_final_condition_spectra.R --results-dir "$TSF_RESULTS_DIR" --cores "$N_WORKERS"

  step "[$NAME] crest genes per characteristic peak, built atomically"
  crest_all() {  # $1 = out dir (last arg, supplied by ensure_dir_atomic)
    local outdir="${@: -1}" sig cond n
    for sig in "$TSF_RESULTS_DIR"/condition_library/condition_signature_*.tsv; do
      [ -f "$sig" ] || continue
      cond=$(basename "$sig" .tsv); cond=${cond#condition_signature_}
      n=$(($(wc -l < "$sig") - 1)); [ "$n" -gt 0 ] || { log "$cond: empty signature, skipped"; continue; }
      Rscript scripts/crest_genes.R --signature "$sig" --datasets "$COHORTS" --condition "$cond" \
        --out "$outdir" || log "crest genes failed for $cond (continuing)"
    done
  }
  ensure_dir_atomic "$TSF_RESULTS_DIR/crest_genes" "$REF_DIG|lib=$(md5 "$TSF_RESULTS_DIR/condition_library.inputs")" --out crest_all

  step "[$NAME] gene-level LOCO baseline"
  # a failure here must not abort the run: the libraries that follow do not
  # depend on it, and the baseline can be re-run alone afterwards
  python3 scripts/run_gene_baseline.py --interim-dir "$TSF_INTERIM_DIR" --datasets "$COHORTS" \
    --target class_id --out "$TSF_RESULTS_DIR/gene_baseline.tsv" \
    || log "[$NAME] gene baseline FAILED (continuing); re-run: python3 scripts/run_gene_baseline.py --interim-dir $TSF_INTERIM_DIR --datasets $COHORTS --target class_id --out $TSF_RESULTS_DIR/gene_baseline.tsv"
  log "[$NAME] done -> $TSF_RESULTS_DIR"
}

# ------------------------------------------------------------- 4. primary
PRIMARY="$GTEX_ID${R3_COHORTS:+,$R3_COHORTS}"
want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }
if ! want primary; then :
elif [ -z "$R3_COHORTS" ]; then
  log "PRIMARY library would be GTEx alone: no cohort of $TSF_GSE is in recount3, so there is no"
  log "recount3-only disease cohort to validate against. The primary library is skipped; the"
  log "sensitivity library (GTEx + GEO) is the only one buildable, and it MUST be reported as such."
else
  run_library primary "$PRIMARY"
fi

# ---------------------------------------------------------- 5. sensitivity
if want sensitivity_geo && [ -n "$GEO_ONLY" ]; then
  run_library sensitivity_geo "$GTEX_ID,$GEO_ONLY"
fi

# ------------------------------------------------------------- 6. combined
[ "$ONLY" = "combined" ] && TSF_RUN_COMBINED=1
if ! want combined; then :
elif [ "$TSF_RUN_COMBINED" = "1" ] && [ -n "$R3_COHORTS" ] && [ -n "$GEO_ONLY" ]; then
  run_library combined "$GTEX_ID,$R3_COHORTS,$GEO_ONLY"
elif [ -n "$R3_COHORTS" ] && [ -n "$GEO_ONLY" ]; then
  log "combined library not built (TSF_RUN_COMBINED=1 enables it). Build it only after comparing"
  log "the LOCO of primary and sensitivity_geo and seeing that technical origin does not dominate."
fi

# ------------------------------------------------------------- 7. summary
step "summary"
for L in primary sensitivity_geo combined; do
  [ -d "$ROOT_RESULTS/$L" ] || continue
  log "$L:"
  log "   inputs digest   $ROOT_RESULTS/$L/inputs_digest.txt"
  log "   profile         $ROOT_INTERIM/$L/reference_profile_${TSF_TISSUE}.tsv (+ _manifest.tsv)"
  log "   library         $ROOT_RESULTS/$L/condition_library/"
  log "   crest genes     $ROOT_RESULTS/$L/crest_genes/"
  log "   spectral LOCO   $ROOT_RESULTS/$L/reference/"
  log "   gene LOCO       $ROOT_RESULTS/$L/gene_baseline.tsv"
done
log "Compare cohort_drop (within-cohort minus out-of-cohort) of the spectral and gene LOCOs, per library."
log "Match a new sample against a library:  ./tsf match <counts.tsv> --results-dir $ROOT_RESULTS/primary"
log "(the library carries its own profile; TSF_REFERENCE_PROFILE is not needed)"
