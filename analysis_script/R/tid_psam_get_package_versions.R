# tid_psam_get_package_versions.R
#
# Prints out the versions of used packages.
#
# Tim Dressler, 21.06.2025


#-------------------------------------Set up------------------------------------
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Set up paths
if (interactive()) {
  # Running inside RStudio
  library(rstudioapi)
  SCRIPTPATH <- dirname(rstudioapi::getSourceEditorContext()$path)
} else {
  # Running via Rscript (called from tid_psam_run_pipeline.py)
  args <- commandArgs(trailingOnly = FALSE)
  script_path <- normalizePath(sub("--file=", "", args[grep("--file=", args)]))
  SCRIPTPATH <- dirname(script_path)
}
if (grepl("psam/analysis_script/R", SCRIPTPATH)) {
  cat("Path OK") 
} else {
  stop("Path not OK")  
}


# Get packages
extract_packages <- function(file) {
  lines <- readLines(file, warn = FALSE)
  libs <- unlist(regmatches(lines, gregexpr("(?<=library\\()[^\\)]+", lines, perl = TRUE)))
  reqs <- unlist(regmatches(lines, gregexpr("(?<=require\\()[^\\)]+", lines, perl = TRUE)))
  colons <- unlist(regmatches(lines, gregexpr("\\b\\w+(?=::)", lines, perl = TRUE)))
  unique(trimws(c(libs, reqs, colons)))
}

r_files <- list.files(path = SCRIPTPATH, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
all_packages <- unique(unlist(lapply(r_files, extract_packages)))

pkg_versions <- sapply(all_packages, function(pkg) {
  tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
})

# Print versions
cat("Packages used in R scripts and their installed versions:\n\n")
print(data.frame(Package = all_packages, Version = pkg_versions, row.names = NULL))
