# tid_psam_check_folder_TD - Checks if the specified folder paths exist. If a path does not exist, the function creates it.
#
# Usage:
#   tid_psam_check_folder_TD('Path1', 'Path2', ...)
#
# Inputs:
#   ... - any number of folder paths as strings
#
# Tim Dressler, 29.05.2025

tid_psam_check_folder_TD <- function(...) {
  folders <- list(...)
  
  for (folder in folders) {
    if (!dir.exists(folder)) {
      cat(sprintf('The folder "%s" does not exist.\n', folder))
      success <- tryCatch({
        dir.create(folder, recursive = TRUE)
        TRUE
      }, error = function(e) {
        FALSE
      })
      
      if (success) {
        cat(sprintf('Folder "%s" was successfully created.\n', folder))
      } else {
        cat(sprintf('Error: Could not create folder "%s".\n', folder))
      }
    } else {
      cat(sprintf('Folder "%s" exists.\n', folder))
    }
  }
}
