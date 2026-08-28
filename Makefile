.PHONY: test check run selfcheck status clean

test:
	Rscript tests/test_labels.R
	Rscript tests/test_spectrum.R

check:
	./tsf check

run:
	./tsf run

selfcheck:
	./tsf selfcheck

status:
	./tsf status

clean:
	rm -rf $(TSF_RESULTS_DIR)/*
