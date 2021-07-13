# Arun Chavan
# Started: 2021-04-01

# Oraganize file names, checksums, etc, for preparing the data for submission to GEO.

# set up ======================================================================
library(tidyverse)
library(glue)

# read filenames and checksums for fastq files ================================
# These files are in the scratch space on farnam temporarily, from where I'll
# upload them to GEO. I rand md5sum on them there, wrote the results to a file,
# which we are reading here.
files <- read_delim("data/geo/update_20210710/checksums.txt", 
                    col_names = c("md5", "filename"), 
                    delim = " ", 
                    trim_ws = TRUE)

# read sample info ============================================================
sample.info <- read_csv("ext/from-alice/sample_info.csv")

# join ========================================================================
# extract AL id
files$AL_id <- str_extract(files$filename, "AL-*\\d+") %>%
  str_replace("-", "") %>% 
  str_replace(pattern = "(AL)(\\d{1})$", replacement = "\\10\\2") # left pad

files <- inner_join(sample.info, files, by = "AL_id")

# new filenames ===============================================================
# Create new filenames that are based on sample info rather than AL id, so that
# they are informative when we upload them to GEO.
files <- files %>% 
  mutate(filename_new = str_replace(filename, 
                                    "(AL[-]*[\\d]+)(.*)", 
                                    glue("{sample_id}__{tissue}_\\2")))

# save file info
files <- files %>% 
  relocate(md5, .after = filename_new) %>% 
  arrange(sample_id)

write_tsv(files, "data/geo/update_20210710/file-info.tsv") 

# write a file info file with only the rows for the new files for AL8 that will
# be uploaded to SRA for updating the record.
files %>% 
  filter(str_detect(filename, "AL8_HCR_S10")) %>% 
  write_tsv("data/geo/update_20210710/file-info_new-files.tsv")

# script for renaming files ===================================================
# based on the above, create a script for renaming the fastq files
str_c(
  "#!/usr/bin/env bash\n\n", 
  "cd /home/arc78/scratch60/covid-placenta/data/00_fastq\n\n",
  str_flatten(glue("mv {files$filename} {files$filename_new}"), 
              collapse = "\n"), 
  "\n\n"
) %>%
  write_lines("code/00_rename-fastq-files.sh")

# end =========================================================================







