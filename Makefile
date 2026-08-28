.PHONY: test ingest clean

test:
	Rscript tests/test_labels.R

ingest:
	Rscript scripts/01_ingest_dataset.R

clean:
	rm -rf data/interim/*
