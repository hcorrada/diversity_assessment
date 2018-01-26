########################## Generate data frame with seq data set characteristics 
# columns
# - seq run 
# - individual  
# - sample id  
# - titration  
# - number of reads  
# - average base quality 

qa_file <- "data/rqcQA_list.rds"

## Extract read count info from qc output and cleanup data frame
qa_file_info <- readRDS(qa_file) %>%
    perFileInformation() %>%
    select(-format, -path) %>%
    mutate(filename = as.character(filename))

qa_file_info <- qa_file_info %>% select(-pair, -group)
if (all(qa_file_info$reads == qa_file_info$total.reads)) {
    qa_file_info$total.reads <- NULL
}

read_count_df <- qa_file_info %>% 
    mutate(sample_id = str_extract(filename, ".*(?=_S)"),
           read_dir = str_extract(filename, "(?<=L001_)..(?=_001)")) %>% 
    select(-filename) %>% 
    left_join(mgtstMetadata) 

ProjectTemplate::cache("seq_char_df", 
                       CODE = {read_count_df}, 
                       depends = mgtstMetadata)


## Data frame with PhiX error Info 
sav_dir <- "../data/sav"
sav_list <- list(
    jhu1 = file.path(sav_dir,"run_18071097_sav"),
    jhu2 = file.path(sav_dir,"run_18527531_sav"),
    nist1 = file.path(sav_dir, "run_170209_sav"),
    nist2 = file.path(sav_dir, "run_170216_sav")
)

error_df <- sav_list %>% map(savR) %>%
    map_df(errorMetrics,.id = "ds") %>%
    mutate(read = if_else(cycle < 301, "R1","R2")) %>% 
    mutate(error_df, base_position = if_else(cycle < 301, cycle, (450 - cycle) + 300)) 

ProjectTemplate::cache("phix_error_df", 
                       CODE = {error_df}, 
                       depends = mgtstMetadata)

## Pipeline characterization
## Columns
## - Number of features 
## - total abundance 
## - filter rate

# cache("seq_and_pipe_info")
