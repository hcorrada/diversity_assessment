## Calculate signal to noise ratio comparing titrations and unmixed post to unmixed pre

############ Within PCR replicate pairwise distance summary

get_within_dist_summary_df <- function(dist_obj, mgtstMetadata){
    ## pairwise index combinations
    idx_combn <- combn(1:4, 2,simplify = FALSE)
    
    ## sample names and number of samples from distance object
    dist_labels <- labels(dist_obj)
    sams <- length(dist_labels)
    
    dist_idx_df <- mgtstMetadata %>% 
        filter(biosample_id != "NTC") %>% 
        ## Get sample idx value for dist object
        mutate(dist_idx = match(sample_id, dist_labels))
    
    ## Create list of replicate ids 
    rep_idx_df <- dist_idx_df %>% 
        group_by(seq_lab, seq_run, biosample_id, t_fctr) %>% 
        select(-pcr_16S_plate, -pos) %>% 
        nest() %>% 
        mutate(rep_idx = map(data, ~.$dist_idx))
    
    ## replicate idx combinations
    rep_pair_idx_df <- rep_idx_df %>% 
        mutate(idx_combn = list(idx_combn)) %>% 
        unnest(idx_combn,.drop = FALSE) %>% 
        mutate(rep_comb = map2(rep_idx, idx_combn, ~.x[.y]))
    
    ## extract pairdist value from dist_object   
    rep_pair_dist_df <-  rep_pair_idx_df  %>% 
        mutate(rep_dist = map_dbl(rep_comb, ~dist_obj[sams*(.[1]-1) - .[1]*(.[1]-1)/2 + .[2]-.[1]]))
    
    ## pairwise distance summary statistics
    rep_pair_dist_df %>% 
        group_by(seq_lab, seq_run, biosample_id, t_fctr) %>% 
        summarise(mean_rep_dist = mean(rep_dist),
                  median_rep_dist = median(rep_dist),
                  max_rep_dist = max(rep_dist),
                  min_rep_dist = min(rep_dist), 
                  sd_rep_dist = sd(rep_dist))
}

############### Between titration pairwise distance summary 
get_between_dist_summary_df <- function(dist_obj, mgtstMetadata){
    ## sample names and number of samples from distance object
    dist_labels <- labels(dist_obj)
    sams <- length(dist_labels)
        
    dist_idx_df <- mgtstMetadata %>% 
        filter(biosample_id != "NTC") %>% 
        ## Get sample idx value for dist object
        mutate(dist_idx = match(sample_id, dist_labels))
    
    pre_dist_idx_df <- dist_idx_df %>% 
        filter(t_fctr == 20) %>% 
        dplyr::rename(pre_idx = dist_idx, pre_id = sample_id) %>% 
        select(-t_fctr, -pcr_16S_plate, -pos)
    
    comp_dist_idx_df <- dist_idx_df %>% filter(t_fctr != 20)
    
    pre_comp_dist_idx_df <- left_join(comp_dist_idx_df, pre_dist_idx_df)
    
    ## Comparison distance summary
    pre_comp_dist_idx_df %>% 
        mutate(pair_dist = dist_obj[sams*(pre_idx - 1) - pre_idx*(pre_idx - 1)/2 + dist_idx - pre_idx]) %>% 
        group_by(seq_lab, seq_run, biosample_id, t_fctr) %>% 
        summarise(mean_dist = mean(pair_dist),
                  median_dist = median(pair_dist),
                  max_dist = max(pair_dist),
                  min_dist = min(pair_dist), 
                  sd_dist = sd(pair_dist))
}

############### Signal to noise df 
get_signal_to_noise_df <- function(within_dist_summary_df, between_dist_summary_df){
    pre_within_dist <- within_dist_summary_df %>% 
        filter(t_fctr == 20) %>% 
        dplyr::rename(mean_pre_dist = mean_rep_dist, 
                      median_pre_dist = median_rep_dist) %>% 
        select(-contains("rep_dist"), -t_fctr)
    
    between_dist_summary_df %>% 
        left_join(within_dist_summary_df) %>% 
        left_join(pre_within_dist) %>% 
        rowwise() %>% 
        mutate(mean_noise = mean(c(mean_rep_dist, mean_pre_dist),na.rm = TRUE), 
               median_noise = mean(c(mean_rep_dist, mean_pre_dist), na.rm = TRUE),
               mean_sig_noise = mean_dist/mean_noise,
               median_sig_noise = median_dist/median_noise)
}

############### Signal to noise analysis 
run_signal_to_noise <- function(dist_methods, mgtstMetadata){  
    ## Get dist data frame  
    beta_df <- make_beta_div_df(dist_methods) %>% 
        mutate(dist_obj = map(dist_results, pluck, "result")) %>%
        filter(!is.null(dist_obj)) %>% 
        select(-dist_results)
    
    if (dist_methods == "jaccard") {
        beta_df <- beta_df %>% 
            mutate(dist_obj = map(dist_obj, ~(1 - .)))
    }
    
    ## within pair distance summary
    within_dist_summary_df <- beta_df %>%
        mutate(pair_df = map(dist_obj, 
                             get_within_dist_summary_df, 
                             mgtstMetadata)) %>% 
        select(-dist_obj) %>%
        unnest()
    
    ## Between pair distance summary   
    between_dist_summary_df <- beta_df %>%
        mutate(pair_df = map(dist_obj, 
                             get_between_dist_summary_df, 
                             mgtstMetadata)) %>% 
        select(-dist_obj) %>%
        unnest()

    ## Signal to noise summary 
    get_signal_to_noise_df(within_dist_summary_df, between_dist_summary_df)
}

######################## Bray Curtis Evaluation ################################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "sig_noise_bray_df",
    {
        run_signal_to_noise("bray", mgtstMetadata)
    },
    depends = c("mgtstMetadata")
)

######################## Weighted UniFrac Evaluation ###########################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "sig_noise_wunifrac_df",
    {
        run_signal_to_noise("wunifrac", mgtstMetadata)
    },
    depends = c("mgtstMetadata")
)

######################## Jaccard Evaluation ####################################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "sig_noise_jaccard_df",
    {
        run_signal_to_noise("jaccard", mgtstMetadata)
    },
    depends = c("mgtstMetadata")
)

######################## Unifrac Metric Evaluation #############################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "sig_noise_unifrac_df", 
    {
        run_signal_to_noise("unifrac", mgtstMetadata)
    }, 
    depends = c("mgtstMetadata")
)
