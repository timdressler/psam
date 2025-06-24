# tid_psam_clean_up_folder_TD - Moves files to a timestamped archive subfolder within an "archive" folder
#
# Usage:
#   tid_psam_clean_up_folder_TD('Path')
#
# Input:
#   folderPath - The path of the folder to check and clean up.
#
# Description:
#   This function checks if a specified folder is empty. If the folder
#   contains any files, it ensures an "archive" folder exists within the
#   specified folder, creates a new subfolder named archive_DD_MM_YYYY_HH-MM-SS
#   inside the archive folder, and moves all files from the original folder
#   to this new archive subfolder.
#
# Tim Dressler, 29.05.2025

tid_psam_clean_up_folder_TD <- function(folder_path) {
  if (!dir.exists(folder_path)) {
    stop("Specified folder does not exist.")
  }
  
  items <- list.files(folder_path, full.names = TRUE, all.files = FALSE, recursive = FALSE)
  
  # Exclude "archive" 
  items <- items[!basename(items) %in% c("archive")]
  
  if (length(items) == 0) {
    cat("The folder is empty. No action taken.\n")
    return()
  }
  
  archive_folder <- file.path(folder_path, "archive")
  if (!dir.exists(archive_folder)) {
    dir.create(archive_folder)
  }
  
  timestamp <- format(Sys.time(), "%d_%m_%Y_%H-%M-%S")
  archive_subfolder <- file.path(archive_folder, paste0("archive_", timestamp))
  dir.create(archive_subfolder)
  
  for (item in items) {
    dest <- file.path(archive_subfolder, basename(item))
    file.rename(item, dest)
  }
  
  cat(sprintf("All contents moved to: %s\n", archive_subfolder))
}
