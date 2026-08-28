.PHONY: test check ingest spectra maxt stability peaks compare all clean

test:
	Rscript tests/test_labels.R
	Rscript tests/test_spectrum.R

check:
	Rscript scripts/00_check_inputs.R

ingest: check
	Rscript scripts/01_ingest_dataset.R

spectra:
	Rscript scripts/02_spectra.R

maxt:
	Rscript scripts/03_maxt.R

stability:
	Rscript scripts/04_stability.R

peaks:
	Rscript scripts/05_peaks_genes.R

compare:
	Rscript scripts/06_compare_datasets.R

all: ingest spectra maxt stability peaks compare

clean:
	rm -rf $(TSF_RESULTS_DIR)/*
