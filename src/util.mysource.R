mysource <- function(filename) {
  filelocations <- c(filename, paste0("src/", filename), paste0("../src/", filename))
  warning_msg_if_not_found <- paste0("Make sure ", filename, " is sourced; otherwise error may occur.")
  
  for (f in filelocations) {
    if (file.exists(f)) {
      source(f)
      return(0)
    }
  }
  stop(warning_msg_if_not_found)
}


