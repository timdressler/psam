# tid_psam_install_requirements_R.R
#
# Installs needed R packages if they are not already installed.
#
# Tim Dressler, 16.06.25

# List of all required packages 
required_packages <- unique(c(
  "tidyr", "afex", "emmeans", "readxl", "car", "corrplot", "dplyr", 
  "ez", "ggplot2", "ggstatsplot", "DescTools", "ggpubr", "ggeffects", 
  "cowplot", "tidyverse", "psych", "rstatix", "lmerTest", "RVAideMemoire" 
))

# Get currently installed packages
installed_packages <- installed.packages()[, "Package"]

# Track which ones are newly installed
newly_installed <- c()

# Install missing packages
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, dependencies = TRUE)
    newly_installed <- c(newly_installed, pkg)
  }
}

# Display the newly installed packages (if any)
if (length(newly_installed) > 0) {
  cat("The following packages were installed:\n")
  print(newly_installed)
} else {
  cat("All required packages were already installed.\n")
}
