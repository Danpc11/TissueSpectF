.PHONY: test test-r test-ml check run selfcheck status clean clean-dry clean-force sonify

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
test-ml:
	@if python3 -c 'import numpy, sklearn, pytest' 2>/dev/null; then \
	  python3 -m pytest tests/ml/ -q; \
	else \
	  echo "skipping tests/ml: numpy, scikit-learn or pytest missing (pip install -r requirements-ml.txt)"; \
	fi

check:
	./tsf check

run:
	./tsf run

selfcheck:
	./tsf selfcheck

status:
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
