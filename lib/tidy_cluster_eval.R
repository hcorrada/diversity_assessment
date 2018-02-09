## Tidy cluster evaluation results - not returning tidy objects for debugging
tidy_cluster_eval <- function(cluster_eval_df) {
    ## Filtering failed evaluations
    eval_df <- cluster_eval_df %>%
        mutate(eval_result = map(eval_output, pluck, "result")) %>%
        mutate(eval_error = map_lgl(eval_result, is.null)) %>%
        filter(!eval_error) %>% 
        select(pipe, method, dist_method, eval_result) 
    
    ## expanding results and getting number of samples per comparison
    eval_df <- eval_df %>%
        unnest() %>% 
        mutate(n_samples = map_int(cluster_output, ~length(.$result))) 
    
    ## Tidying results %>%
    eval_df %>%  
        select(-comp_df, -cluster_output) %>% 
        ## Excluding NULL results - comparisons with less than 2 samples 
        filter(cluster_results != -1)
}
