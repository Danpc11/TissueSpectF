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

# Deletes the results tree named in config/project.R -- never an unguarded path.
clean:
	@Rscript -e 'source("R/utils_io.R"); source("R/config.R"); \
	  p <- load_project_config("config/project.R"); d <- p$$results_dir; \
	  if (is.null(d) || !nzchar(d) || d %in% c("/", "/*") || nchar(d) < 8) \
	    stop("Refusing to clean unsafe path: ", d, call. = FALSE); \
	  if (!dir.exists(d)) { cat("nothing to clean:", d, "\n"); quit() }; \
	  cat("removing contents of", d, "\n"); \
	  unlink(list.files(d, full.names = TRUE), recursive = TRUE)'
