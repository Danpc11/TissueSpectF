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
# everything that went into it (dataset list, GTEx project, tissue, GSE list,
# annotation, gene universe, axis, bin size, expression filters, profile
# digest). The artefact is rebuilt when that digest differs, so changing
# TSF_GSE or TSF_BIN_SIZE cannot silently reuse a mask or a profile built for
# another configuration. Re-run the same command after a failure and it
# continues.
#
# Required environment:
#   TSF_ROOT  TSF_GEO_DIR  TSF_INTERIM_DIR  TSF_RESULTS_DIR   (interim/results: NEW trees)
# Optional:
#   TSF_GSE          GEO cohorts (default: the five liver ones)
#   TSF_GTEX_TISSUE  recount3 GTEx project (LIVER)      TSF_TISSUE  label (liver)
#   TSF_VOCAB        vocabulary (liver_fibrosis; must contain Control_external_study)
#   TSF_BIN_SIZE     bp bin width (100000)              N_WORKERS   cores (4)
#   TSF_RUN_COMBINED=1  also build the combined library
#   TSF_CONDITION_B / TSF_MAXT_B  permutations (defaults from config/project.R)
#
# Usage:  bash scripts/run_differential.sh 2>&1 | tee "$TSF_RESULTS_DIR/run_differential.log"
set -euo pipefail

: "${TSF_ROOT:?export TSF_ROOT}"; : "${TSF_GEO_DIR:?export TSF_GEO_DIR}"
: "${TSF_INTERIM_DIR:?export TSF_INTERIM_DIR}"; : "${TSF_RESULTS_DIR:?export TSF_RESULTS_DIR}"
TSF_GSE="${TSF_GSE:-GSE135251,GSE130970,GSE162694,GSE276114,GSE142530}"
TSF_GTEX_TISSUE="${TSF_GTEX_TISSUE:-LIVER}"
TSF_TISSUE="${TSF_TISSUE:-liver}"
TSF_VOCAB="${TSF_VOCAB:-liver_fibrosis}"
TSF_BIN_SIZE="${TSF_BIN_SIZE:-100000}"
N_WORKERS="${N_WORKERS:-4}"
TSF_RUN_COMBINED="${TSF_RUN_COMBINED:-0}"
cd "$TSF_ROOT"
ROOT_INTERIM="$TSF_INTERIM_DIR"; ROOT_RESULTS="$TSF_RESULTS_DIR"
mkdir -p "$TSF_GEO_DIR" "$ROOT_INTERIM" "$ROOT_RESULTS"
log()  { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
step() { log "=== $* ==="; }
md5()  { md5sum "$1" | cut -c1-32; }

# ---------------------------------------------------------------- inputs digest
# The configuration every artefact depends on. Anything here changes -> rebuild.
ANNOT=$(grep -o 'annotation_file *= *[^,]*' config/project.R | head -1 | tr -d ' ')
UNIVERSE="${TSF_GENE_UNIVERSE:-$(grep -o 'gene_universe *= *[^,]*' config/project.R | head -1 | tr -d ' ')}"
FILTERS=$(grep -oE '(min_tpm|min_fraction) *= *[0-9.]+' config/project.R | tr -d ' ' | paste -sd, -)
PROJECT_MD5=$(md5 config/project.R)
base_digest() {  # $1 = dataset list
  printf 'datasets=%s|gtex=%s|tissue=%s|gse=%s|annot=%s|universe=%s|axis=bp|bin=%s|filters=%s|project=%s' \
    "$1" "$TSF_GTEX_TISSUE" "$TSF_TISSUE" "$TSF_GSE" "$ANNOT" "$UNIVERSE" "$TSF_BIN_SIZE" "$FILTERS" "$PROJECT_MD5"
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

# ---------------------------------------------------------------- 0. sanity
step "0 tests and recount3 config validation"
TESTS_DIG="tests=$(git rev-parse HEAD 2>/dev/null || echo nogit)"
ensure "$ROOT_RESULTS/.tests_ok" "$TESTS_DIG" bash -c 'make test >/dev/null && touch "$0"' "$ROOT_RESULTS/.tests_ok"
if ls config/datasets/R3_*.R >/dev/null 2>&1; then
  Rscript scripts/validate_recount3_configs.R --vocabulary "$TSF_VOCAB" || {
    log "STOP: recount3 configs from an earlier generator. Review them, then:"
    log "      Rscript scripts/validate_recount3_configs.R --migrate --vocabulary $TSF_VOCAB"
    log "      (header fields only; condition_rules are never touched; .bak kept)"
    exit 2
  }
fi

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
ensure "$TSF_GEO_DIR/${GTEX_ID}_reads.tsv.gz" "gtex=$TSF_GTEX_TISSUE|tissue=$TSF_TISSUE|vocab=$TSF_VOCAB" \
  Rscript scripts/recount3_fetch.R --projects "$TSF_GTEX_TISSUE" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
R3_COHORTS=""
if [ -n "$SRPS" ]; then
  for s in ${SRPS//,/ }; do
    ensure "$TSF_GEO_DIR/R3_${s}_reads.tsv.gz" "srp=$s|tissue=$TSF_TISSUE|vocab=$TSF_VOCAB" \
      Rscript scripts/recount3_fetch.R --projects "$s" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
    if grep -q '^\s*# list(id = "biopsy_fibrosis_stage"' "config/datasets/R3_$s.R"; then
      log "STOP: config/datasets/R3_$s.R needs condition_rules for the exploded SRA attributes (the pheno columns are listed in it). Edit it, then re-run."
      exit 2
    fi
    R3_COHORTS="${R3_COHORTS:+$R3_COHORTS,}R3_$s"
  done
  Rscript scripts/validate_recount3_configs.R --vocabulary "$TSF_VOCAB" || exit 2
fi

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
  ensure "$TSF_RESULTS_DIR/.stages_done" "$DIG|mask=$MASK_MD5|profile=$PROF_MD5|B=${TSF_CONDITION_B:-cfg}/${TSF_MAXT_B:-cfg}" \
    bash -c './tsf run $1 --from spectra --to compare --grid-axis bp --bin-size "$2" --gene-mask "$3" --stage-order F0,F1,F2,F3,F4 && touch "$4"' _ \
      "$DS_LIST" "$TSF_BIN_SIZE" "$MASK" "$TSF_RESULTS_DIR/.stages_done"
  step "[$NAME] fingerprint library + out-of-cohort validation (profile stored in reference.rds)"
  ensure "$TSF_RESULTS_DIR/reference/reference.rds" "$DIG|mask=$MASK_MD5|profile=$PROF_MD5" \
    ./tsf reference $DS_LIST --grid-axis bp --bin-size "$TSF_BIN_SIZE" --gene-mask "$MASK"

  step "[$NAME] condition library (cross-cohort meta-analysis)"
  ensure "$TSF_RESULTS_DIR/condition_library" "$DIG|mask=$MASK_MD5|profile=$PROF_MD5" \
    Rscript scripts/build_final_condition_spectra.R --results-dir "$TSF_RESULTS_DIR" --cores "$N_WORKERS" \
      --out-dir "$TSF_RESULTS_DIR/condition_library"

  step "[$NAME] crest genes per characteristic peak"
  for sig in "$TSF_RESULTS_DIR"/condition_library/condition_signature_*.tsv; do
    [ -f "$sig" ] || continue
    local cond; cond=$(basename "$sig" .tsv); cond=${cond#condition_signature_}
    local n; n=$(($(wc -l < "$sig") - 1)); [ "$n" -gt 0 ] || { log "$cond: empty signature, skipped"; continue; }
    Rscript scripts/crest_genes.R --signature "$sig" --datasets "$COHORTS" --condition "$cond" \
      --out "$TSF_RESULTS_DIR/crest_genes" || log "crest genes failed for $cond (continuing)"
  done

  step "[$NAME] gene-level LOCO baseline"
  python3 scripts/run_gene_baseline.py --interim-dir "$TSF_INTERIM_DIR" --datasets "$COHORTS" \
    --target class_id --out "$TSF_RESULTS_DIR/gene_baseline.tsv"
  log "[$NAME] done -> $TSF_RESULTS_DIR"
}

# ------------------------------------------------------------- 4. primary
PRIMARY="$GTEX_ID${R3_COHORTS:+,$R3_COHORTS}"
if [ -z "$R3_COHORTS" ]; then
  log "PRIMARY library would be GTEx alone: no cohort of $TSF_GSE is in recount3, so there is no"
  log "recount3-only disease cohort to validate against. The primary library is skipped; the"
  log "sensitivity library (GTEx + GEO) is the only one buildable, and it MUST be reported as such."
else
  run_library primary "$PRIMARY"
fi

# ---------------------------------------------------------- 5. sensitivity
if [ -n "$GEO_ONLY" ]; then
  run_library sensitivity_geo "$GTEX_ID,$GEO_ONLY"
fi

# ------------------------------------------------------------- 6. combined
if [ "$TSF_RUN_COMBINED" = "1" ] && [ -n "$R3_COHORTS" ] && [ -n "$GEO_ONLY" ]; then
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
log "Match a new sample against a library:  TSF_RESULTS_DIR=$ROOT_RESULTS/primary ./tsf match <counts.tsv>"
log "(the library carries its own profile; TSF_REFERENCE_PROFILE is not needed)"
