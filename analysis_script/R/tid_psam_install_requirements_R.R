# tid_psam_install_requirements_R.R
#
# Installs needed R packages if not already installed.
# Also checks that installed versions match the required ones.
#
# Tim Dressler, 16.06.2025

# Required packages and their used versions
required_packages <- list(
  tidyr = "1.3.1",
  afex = "1.4.1",
  emmeans = "1.11.1",
  readxl = "1.4.5",
  car = "3.1.3",
  corrplot = "0.95",
  dplyr = "1.1.4",
  ez = "4.4.0",
  ggplot2 = "3.5.2",
  ggstatsplot = "0.13.1",
  DescTools = "0.99.60",
  ggpubr = "0.6.0",
  ggeffects = "2.2.1",
  cowplot = "1.1.3",
  tidyverse = "2.0.0",
  psych = "2.5.3",
  rstatix = "0.7.2",
  lmerTest = "3.1.3",
  RVAideMemoire = "0.9.83.11",
  effectsize = "1.0.1",
  smplot2 = "0.2.5",
  nparLD = "2.2",
  Rmisc = "1.5.1",
  devtools = "2.4.5",
  renv = "1.1.4"
)

# Get installed packages and versions
installed <- installed.packages()

# Track newly installed and version mismatches
newly_installed <- c()
version_mismatches <- c()

for (pkg in names(required_packages)) {
  required_version <- required_packages[[pkg]]
  if (!(pkg %in% rownames(installed))) {
    install.packages(pkg, dependencies = TRUE)
    newly_installed <- c(newly_installed, paste0(pkg, " (", required_version, ")"))
  } else {
    current_version <- as.character(packageVersion(pkg))
    if (current_version != required_version) {
      version_mismatches <- c(version_mismatches,
                              paste0(pkg, " (installed: ", current_version, 
                                     ", required: ", required_version, ")"))
    }
  }
}

if (length(newly_installed) > 0) {
  cat("The following packages were newly installed:\n")
  print(newly_installed)
} else {
  cat("All required packages were already installed.\n")
}

if (length(version_mismatches) > 0) {
  cat("\n The following packages are installed but do NOT match the required version:\n")
  print(version_mismatches)
} else {
  cat("\n All installed packages match the required versions.\n")
}
