library(tidyverse)

demopath <- "demopath"
demofile <- file.path(demopath, "DemodfScreenFinal.csv")
cbclfile <- file.path(demopath, "mh_p_cbcl.csv")

demodf <- read.csv(demofile)
cbcl <- read.csv(cbclfile) %>%
  select(src_subject_id, eventname, cbcl_scr_syn_totprob_r)

merged <- demodf %>%
  left_join(cbcl, by = c("src_subject_id", "eventname"))

message("Rows before merge: ", nrow(demodf))
message("Rows after merge:  ", nrow(merged))
message("CBCL total raw NA rate: ", mean(is.na(merged$cbcl_scr_syn_totprob_r)))

write.csv(merged, demofile, row.names = FALSE)
