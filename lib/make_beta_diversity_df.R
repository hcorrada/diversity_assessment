## Generating Data Frame with weighted and unweighted diversity metrics 

#' Make Tidy Beta Diversity Data Frame
#'
#' @param dist_methods character vector with diversity metrics to include in the data frame
#' @param div_rds_path path to RDS files with dist objects for diversity metrics 
#'
#' @return data frame with a list column of diversity metrics
#' @export
#'
#' @examples
make_beta_div_df <- function(dist_methods, div_rds_path = "data/diversity_data"){
    diversity_rds <- list.files(div_rds_path, full.names = TRUE)
    
    diversity_names <- basename(diversity_rds) %>% str_replace(".rds","") 
    
    diversity_df <- data_frame(div_info = diversity_names, 
                               rds_file = diversity_rds) %>% 
        separate(div_info, c("pipe","method","dist_method")) %>% 
        mutate(dist_method = paste0(dist_method, "_dist"))
    
    dist_methods <- paste0(dist_methods, "_dist") 
    
    diversity_df %>%
            filter(dist_method %in% dist_methods) %>%
            mutate(dist_results = map(rds_file, readRDS)) %>% 
            select(-rds_file)
}
