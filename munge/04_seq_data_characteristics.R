########################## Generate data frame with seq data set characteristics 
# columns 
# - `sample_id` - unique id for each sequenced sample   
# - `reads`` - number of reads  
# - `read_dir` - read direction   
# - `quality` - mode quality score  
# - `density` - proportion of reads with mode quality score 
get_seq_info <- function(qa_file){

    qa_list <- readRDS(qa_file)

    read_count_df <- qa_list[['readCounts']] %>% 
        select(-filter, -aligned) %>% 
        rownames_to_column(var = "lane")

    read_qual_df <- qa_list[['readQualityScore']] %>% 
        group_by(lane) %>% top_n(1, density) %>% 
        select(-type) %>% 
        ungroup()
    
    ## data frame with read count and quality info 
    left_join(read_count_df, read_qual_df) %>% 
        mutate(sample_id = str_extract(lane, ".*(?=_S)"),
               read_dir = str_extract(lane, "(?<=L001_)..(?=_001)")) %>% 
        select(-lane)
}

ProjectTemplate::cache("seq_char_df", CODE = {get_seq_info("data/qa_list.rds")})


## Data frame with PhiX error Info 
get_error_df <- function() {
    require(savR)
    sav_dir <- "data/sav"
    sav_list <- list(
        jhu1 = file.path(sav_dir, "run_18071097_sav"),
        jhu2 = file.path(sav_dir, "run_18527531_sav"),
        nist1 = file.path(sav_dir, "run_170209_sav"),
        nist2 = file.path(sav_dir, "run_170216_sav")
    )
    
    sav_list %>% map(savR) %>%
        map_df(errorMetrics, .id = "ds") %>%
        mutate(read = if_else(cycle < 301, "R1", "R2")) %>%
        mutate(base_position = if_else(cycle < 301, cycle, (450 - cycle) + 300))
}

ProjectTemplate::cache("phix_error_df", CODE = {get_error_df()}, depends = "mgtstMetadata")