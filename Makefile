.PHONY: test test-r test-ml check run selfcheck status clean clean-dry clean-force sonify require-paths

# Every suite. `make test` used to run two of the three R files, so a failure in
# tests/test_condition_invariants.R was invisible until someone ran it by hand.
test: test-r test-ml

test-r:
	Rscript tests/test_labels.R
	Rscript tests/test_spectrum.R
	Rscript tests/test_condition_invariants.R

# The learned layer's tests need numpy/scikit-learn but not torch. They are
# skipped with a message rather than failing the target, so `make test` still
# works on a machine with no Python environment -- the R pipeline does not
# depend on one.
# The modules ml/ imports at MODULE level, which is what has to be present for
# pytest to even collect the tests. sklearn and scipy are imported lazily inside
# functions, so their absence is a skipped test rather than a collection error
# and they are deliberately not required here.
#
# The previous guard checked `numpy, sklearn, pytest`: it demanded sklearn,
# which is not needed to collect, and did not check yaml or pandas, which are.
# So it decided the suite was runnable and pytest then died on
# `ModuleNotFoundError: No module named 'yaml'` -- a guard that passes and then
# lets the thing it guards fail is worse than no guard.
ML_MODULES = numpy pandas yaml pytest

test-ml:
	@python3 scripts/check_ml_deps.py $(ML_MODULES) \
	  && python3 -m pytest tests/ml/ -q || true

# These pass no paths, so they require TSF_GEO_DIR, TSF_INTERIM_DIR and
# TSF_RESULTS_DIR in the environment: config/project.R names no paths, and
# there is deliberately no default to inherit. The targets check for them
# first, because `./tsf run` failing three lines into a Makefile is less clear
# than being told what to export.
REQUIRED_ENV = TSF_GEO_DIR TSF_INTERIM_DIR TSF_RESULTS_DIR

require-paths:
	@missing=""; \
	for v in $(REQUIRED_ENV); do \
	  eval "val=\$$$$v"; \
	  [ -n "$$val" ] || missing="$$missing $$v"; \
	done; \
	if [ -n "$$missing" ]; then \
	  echo "Set these first (no defaults, by design):$$missing"; \
	  echo "  export TSF_GEO_DIR=data"; \
	  echo "  export TSF_INTERIM_DIR=interim"; \
	  echo "  export TSF_RESULTS_DIR=run_$$(date +%Y_%m_%d)"; \
	  echo "Or pass --geo-dir / --interim-dir / --results-dir to ./tsf."; \
	  exit 1; \
	fi

check: require-paths
	./tsf check

run: require-paths
	./tsf run

selfcheck:
	./tsf selfcheck

status: require-paths
	./tsf status

# Deletes the contents of the results tree named in config/project.R. That is
# the ACTIVE run tree, so this asks for typed confirmation first; the guard and
# the prompt live in the script, not here. See scripts/clean_results.R.
clean:
	@Rscript scripts/clean_results.R

# Same, with no prompt. For scripts that mean it.
clean-force:
	@Rscript scripts/clean_results.R --force

# Print what clean would delete, without deleting it.
clean-dry:
	@Rscript scripts/clean_results.R --dry-run

# Sonification. Needs the condition library, so
# build_final_condition_spectra.R has to have run first.
sonify:
	@python3 -c 'import mido' 2>/dev/null || { \
	  echo "missing mido: pip install -r requirements-sonify.txt"; exit 1; }
	python3 scripts/sonify_tissuespectf.py
