#!/usr/bin/env bash
# run_differential.sh -- the differential pipeline, end to end.
#
#   GTEx liver (recount3)  ->  tissue reference profile
#   cohorts (recount3 when available, GEO otherwise)  ->  deviation from it
#   spectra / condition / consensus / library on the deviation, bp axis
#   crest genes per characteristic peak
#   gene-level LOCO baseline next to the spectral LOCO
#
# Everything is resumable: each step is skipped when its output exists, and
# ./tsf run keeps its own per-stage outputs. Re-run the same command after a
# failure and it continues.
#
# Required environment (export or put in a .env you source first):
#   TSF_ROOT          repo checkout
#   TSF_GEO_DIR       raw inputs (GEO downloads + recount3 tables)
#   TSF_INTERIM_DIR   the common format; USE A NEW ONE for the differential run
#   TSF_RESULTS_DIR   results;          USE A NEW ONE as well
# Optional:
#   TSF_GSE           GEO cohorts to resolve/use   (default: the five liver ones)
#   TSF_GTEX_TISSUE   recount3 GTEx project        (default: LIVER)
#   TSF_TISSUE        tissue label                  (default: liver)
#   TSF_BIN_SIZE      bp bin width                  (default: 100000)
#   TSF_CONDITION_B / TSF_MAXT_B   permutations (defaults from config/project.R)
#   N_WORKERS         cores for the library builder
#
# Usage:  bash scripts/run_differential.sh 2>&1 | tee "$TSF_RESULTS_DIR/run_differential.log"
set -euo pipefail

: "${TSF_ROOT:?export TSF_ROOT}"; : "${TSF_GEO_DIR:?export TSF_GEO_DIR}"
: "${TSF_INTERIM_DIR:?export TSF_INTERIM_DIR}"; : "${TSF_RESULTS_DIR:?export TSF_RESULTS_DIR}"
TSF_GSE="${TSF_GSE:-GSE135251,GSE130970,GSE162694,GSE276114,GSE142530}"
TSF_GTEX_TISSUE="${TSF_GTEX_TISSUE:-LIVER}"
TSF_TISSUE="${TSF_TISSUE:-liver}"
TSF_VOCAB="${TSF_VOCAB:-liver_fibrosis}"   # must contain Control_external_study (the GTEx label)
TSF_BIN_SIZE="${TSF_BIN_SIZE:-100000}"
N_WORKERS="${N_WORKERS:-4}"
cd "$TSF_ROOT"
mkdir -p "$TSF_GEO_DIR" "$TSF_INTERIM_DIR" "$TSF_RESULTS_DIR"
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
step() { log "=== $* ==="; }

GTEX_ID="R3_${TSF_GTEX_TISSUE}"
PROFILE="$TSF_INTERIM_DIR/reference_profile_${TSF_TISSUE}.tsv"
SOURCES="config/recount3_sources.tsv"
MASK="$TSF_INTERIM_DIR/shared_gene_mask.tsv"

# ---------------------------------------------------------------- 0. sanity
step "0 tests"
[ -f "$TSF_RESULTS_DIR/.tests_ok" ] || { make test >/dev/null && touch "$TSF_RESULTS_DIR/.tests_ok"; }

# ---------------------------------------------------------- 1. which sources
step "1 resolve GEO -> SRP -> recount3"
[ -f "$SOURCES" ] || Rscript scripts/resolve_recount3.R --gse "$TSF_GSE" --out "$SOURCES"
# SRPs recount3 holds:
SRPS=$(awk -F'\t' 'NR>1 && $3=="TRUE" {print $2}' "$SOURCES" | paste -sd, -)
# GSEs that stay on the GEO path:
GEO_ONLY=$(awk -F'\t' 'NR>1 && $3!="TRUE" {print $1}' "$SOURCES" | paste -sd, -)
log "recount3 cohorts: ${SRPS:-none}   GEO-only cohorts: ${GEO_ONLY:-none}"
[ -z "$GEO_ONLY" ] || log "NOTE: GEO-only cohorts are quantified with a different pipeline than the GTEx reference; their deviation includes that difference. Say so in the write-up."

# ------------------------------------------------------------- 2. recount3
step "2 fetch recount3 ($GTEX_ID${SRPS:+, $SRPS})"
[ -f "$TSF_GEO_DIR/${GTEX_ID}_reads.tsv.gz" ] || \
  Rscript scripts/recount3_fetch.R --projects "$TSF_GTEX_TISSUE" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
if [ -n "$SRPS" ]; then
  Rscript scripts/recount3_fetch.R --projects "$SRPS" --tissue "$TSF_TISSUE" --vocabulary "$TSF_VOCAB"
  for s in ${SRPS//,/ }; do
    if grep -q '^\s*# list(id = "biopsy_fibrosis_stage"' "config/datasets/R3_$s.R"; then
      log "STOP: config/datasets/R3_$s.R needs condition_rules for the exploded SRA attributes (see the pheno columns listed in it). Edit it, then re-run."
      exit 2
    fi
  done
fi

# ---------------------------------------------------------- 3. GEO inputs
if [ -n "$GEO_ONLY" ]; then
  step "3 GEO fetch/check for $GEO_ONLY"
  ./tsf fetch ${GEO_ONLY//,/ } || true
  ./tsf check ${GEO_ONLY//,/ }
fi

# ------------------------------------------------------------- 4. ingest
COHORTS="$GTEX_ID"
[ -z "$SRPS" ]     || COHORTS="$COHORTS,$(echo "$SRPS" | sed 's/\([^,]*\)/R3_\1/g')"
[ -z "$GEO_ONLY" ] || COHORTS="$COHORTS,$GEO_ONLY"
log "datasets: $COHORTS"
DS_LIST=${COHORTS//,/ }

step "4a ingest, pass 1 (per-cohort expression filter) -> retained_genes"
for d in $DS_LIST; do
  [ -f "$TSF_INTERIM_DIR/$d/retained_genes.tsv" ] || ./tsf ingest "$d" --grid-axis bp --bin-size "$TSF_BIN_SIZE"
done
step "4b shared gene mask across all cohorts"
[ -f "$MASK" ] || Rscript scripts/shared_gene_mask.R --interim-dir "$TSF_INTERIM_DIR" --datasets "$COHORTS" --out "$MASK"
step "4c ingest, pass 2 with the shared mask"
MASK_MD5=$(md5sum "$MASK" | cut -c1-32)
STAMP="mask=$MASK_MD5 bin=$TSF_BIN_SIZE annot=$(grep -o 'annotation_file *= *"[^"]*"' config/project.R | head -1)"
for d in $DS_LIST; do
  # the marker records WHAT was used, not just that something was: a new mask,
  # bin size or annotation invalidates it
  if [ ! -f "$TSF_INTERIM_DIR/$d/.masked" ] || [ "$(cat "$TSF_INTERIM_DIR/$d/.masked")" != "$STAMP" ]; then
    ./tsf ingest "$d" --grid-axis bp --bin-size "$TSF_BIN_SIZE" --gene-mask "$MASK" --force
    printf '%s' "$STAMP" > "$TSF_INTERIM_DIR/$d/.masked"
  fi
done

# ------------------------------------------------- 5. reference & deviation
step "5 tissue reference profile from $GTEX_ID"
[ -f "$PROFILE" ] || Rscript scripts/build_tissue_reference.R --datasets "$GTEX_ID" --tissue "$TSF_TISSUE" --out "$PROFILE"
step "5b every cohort -> deviation from the profile (GTEx included: its deviations are the healthy class)"
# idempotent and stale-safe: re-run after any ingest --force
Rscript scripts/apply_reference_profile.R --datasets "$COHORTS" --profile "$PROFILE"

# ------------------------------------------------------------- 6. spectra
step "6 spectral pipeline on the deviation (null = all, bp axis): spectra ... compare"
./tsf run $DS_LIST --from spectra --to compare --grid-axis bp --bin-size "$TSF_BIN_SIZE" --gene-mask "$MASK"
step "6b fingerprint library + out-of-cohort validation (stores the profile in reference.rds)"
./tsf reference $DS_LIST --grid-axis bp --bin-size "$TSF_BIN_SIZE" --gene-mask "$MASK"

# ------------------------------------------------------ 7. condition library
step "7 condition library (cross-cohort meta-analysis)"
Rscript scripts/build_final_condition_spectra.R --results-dir "$TSF_RESULTS_DIR" --cores "$N_WORKERS" \
  --out-dir "$TSF_RESULTS_DIR/condition_library"

# ---------------------------------------------------------- 8. crest genes
step "8 crest genes per characteristic peak"
for sig in "$TSF_RESULTS_DIR"/condition_library/condition_signature_*.tsv; do
  cond=$(basename "$sig" .tsv); cond=${cond#condition_signature_}
  n=$(($(wc -l < "$sig") - 1)); [ "$n" -gt 0 ] || { log "$cond: empty signature, skipped"; continue; }
  Rscript scripts/crest_genes.R --signature "$sig" --datasets "$COHORTS" --condition "$cond" \
    --out "$TSF_RESULTS_DIR/crest_genes" || log "crest genes failed for $cond (continuing)"
done

# ---------------------------------------------------------- 9. gene baseline
step "9 gene-level LOCO baseline"
python3 scripts/run_gene_baseline.py --interim-dir "$TSF_INTERIM_DIR" --datasets "$COHORTS" \
  --target class_id --out "$TSF_RESULTS_DIR/gene_baseline.tsv"

# ------------------------------------------------------------- 10. summary
step "10 summary"
log "profile:          $PROFILE"
log "library:          $TSF_RESULTS_DIR/condition_library/"
log "crest genes:      $TSF_RESULTS_DIR/crest_genes/"
log "spectral LOCO:    $TSF_RESULTS_DIR/reference/   (validation tables)"
log "gene LOCO:        $TSF_RESULTS_DIR/gene_baseline.tsv"
log "Compare cohort_drop (within-cohort minus out-of-cohort) of the two LOCOs: that is the claim."
log "To match a new sample: ./tsf match <counts.tsv>   (the library carries the profile; no variable needed)"
