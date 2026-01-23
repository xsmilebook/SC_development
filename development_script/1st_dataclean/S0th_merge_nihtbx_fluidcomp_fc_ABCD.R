library(tidyverse)

demopath <- "demopath"
demofile <- file.path(demopath, "DemodfScreenFinal.csv")
nihtbfile <- file.path(demopath, "nc_y_nihtb.csv")

if (!file.exists(demofile)) stop("Missing: ", demofile)
if (!file.exists(nihtbfile)) stop("Missing: ", nihtbfile)

demodf <- read.csv(demofile, stringsAsFactors = FALSE)
nihtb <- read.csv(nihtbfile, stringsAsFactors = FALSE) %>%
  select(src_subject_id, eventname, nihtbx_fluidcomp_fc) %>%
  arrange(src_subject_id, eventname) %>%
  distinct(src_subject_id, eventname, .keep_all = TRUE)

merged <- demodf %>%
  left_join(nihtb, by = c("src_subject_id", "eventname"))

message("Rows before merge: ", nrow(demodf))
message("Rows after merge:  ", nrow(merged))
message("NIHTBX fluidcomp fully-corrected (fc) NA rate: ", mean(is.na(merged$nihtbx_fluidcomp_fc)))

bakfile <- paste0(demofile, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S"))
file.copy(demofile, bakfile, overwrite = FALSE)
message("Backup written: ", bakfile)

write.csv(merged, demofile, row.names = FALSE)
