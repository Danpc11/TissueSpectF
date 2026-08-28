.PHONY: test check ingest clean

test:
	Rscript tests/test_labels.R

check:
	Rscript scripts/00_check_inputs.R

ingest: check
	Rscript scripts/01_ingest_dataset.R

clean:
	rm -rf $(TSF_INTERIM_DIR)/*
