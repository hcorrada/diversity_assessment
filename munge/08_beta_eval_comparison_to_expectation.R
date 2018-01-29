## PAM base cluster evaluation similar to normalization method assessment for
## beta diversity metrics used in McMurdie and Holmes 2014 and Weiss et al 2017
## comparing unmixed pre-exposure samples to titration and post-exposure samples

##################### Functions for PAM Cluster Evaluation #####################

#' Bicluster samples using PAM
#'
#' @param comp_df data frame with samples in sample_ids column
#' @param dist_mat matrix with pairwise distances
#'
#' @return named vector with cluster assignments
#' @export
#'
#' @examples
cluster_samples <- function(comp_df, dist_mat) {
    comp_mat <- dist_mat[row.names(dist_mat) %in% comp_df$sample_id,
                         colnames(dist_mat) %in% comp_df$sample_id]
    if (nrow(comp_mat) <= 2) {
        return(NULL)
    }
    cluster::pam(comp_mat, k = 2, cluster.only = TRUE)
}

#' Evaluate PAM clustering results
#'
#' @param cluster_assignments named vector with cluster assignments
#' @param comp_df data frame with sample_ids (samples) and t_fctr column (grouping)
#'
#' @return double with percent of matched samples
#' @export
#'
#' @examples
cluster_eval <- function(cluster_assignments, comp_df) {
    if (is.null(cluster_assignments)) {
        return(NA)
    }
    
    if (nrow(comp_df) != length(cluster_assignments)) {
        comp_df <-
            filter(comp_df, sample_id %in% names(cluster_assignments))
    }
    
    ## Ensuring sample_ids are in the correct order
    results_df <-
        data.frame(cluster_results = cluster_assignments) %>%
        rownames_to_column(var = "sample_id") %>%
        left_join(comp_df)
    
    ## Number of samples clustered
    samples_compared <- length(cluster_assignments)
    
    ## Evaluate comparison for either possible cluster assignment
    cluster_eval <-
        ifelse(results_df$t_fctr == min(results_df$t_fctr), 1, 2)
    cluster_eval_inv <-
        ifelse(results_df$t_fctr == max(results_df$t_fctr), 1, 2)
    
    assignment_matches <-
        max(
            sum(cluster_eval == results_df$cluster_results),
            sum(cluster_eval_inv == results_df$cluster_results)
        )
    
    ## Return the fraction of samples correctly clustered
    assignment_matches / samples_compared
}

#' Perform PAM cluster evaluation
#'
#' @param dist_obj dist-class object with pairwise distances
#' @param full_comp_df data frame with full set of comparisons
#'
#' @return data frame with cluster evaluation results
#' @export
#'
#' @examples
get_cluster_eval_df <- function(dist_obj, full_comp_df) {
    dist_mat <- as.matrix(dist_obj)
    
    full_comp_df %>%
        mutate(cluster_assignments = map(comp_df, cluster_samples, dist_mat)) %>%
        mutate(cluster_results = map2_dbl(cluster_assignments, comp_df, cluster_eval))
}

## Safe eval version - need to check error after updating source data
safe_get_cluster_eval_df <- safely(get_cluster_eval_df)


perform_cluster_eval <- function(dist_methods, full_comp_df){  
    ## Get dist data frame  
    beta_df <- make_beta_div_df(dist_methods)
    
    ## Perform evaluation
    cluster_eval_df <- beta_df %>%
        mutate(dist_obj = map(dist_results, pluck, "result")) %>%
        filter(!is.null(dist_obj)) %>%
        mutate(eval_output = map(dist_obj, safe_get_cluster_eval_df, full_comp_df))
    
    ## Tidy cluster evaluation results
    cluster_eval_df %>%
        mutate(eval_result = map(eval_output, pluck, "result")) %>%
        mutate(eval_error = map_lgl(eval_result, is.null)) %>%
        filter(!eval_error) %>%
        select(pipe, method, dist_method, eval_result) %>%
        unnest() %>%
        select(-comp_df, -cluster_assignments)
}

######################## Generating Comparison Data Frame ######################
condensed_meta <- mgtstMetadata %>%
    filter(biosample_id != "NTC") %>%
    mutate(sample_id = as.character(sample_id)) %>%
    group_by(seq_lab, seq_run, biosample_id) %>%
    select(-pcr_16S_plate, -pos) %>%
    nest()

comparison_meta <- condensed_meta %>%
    select(-data) %>%
    add_column(t_comp = list(c(0:5, 10, 15))) %>%
    unnest() %>%
    left_join(condensed_meta)

full_comp_df <- comparison_meta %>%
    ## Generating data frame with sample_id and t_fctr for comparisons to perform
    mutate(comp_df = map2(data, t_comp, ~ filter(.x, t_fctr %in% c(.y, 20)))) %>%
    select(-data)

######################## Bray Curtis Evaluation ################################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "beta_cluster_eval_bray_df",
    {
        perform_cluster_eval("bray", full_comp_df)
    },
    depends = c("mgtstMetadata")
)

######################## Weighted UniFrac Evaluation ###########################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "beta_cluster_eval_wunifrac_df",
    {
        perform_cluster_eval("wunifrac", full_comp_df)
    },
    depends = c("mgtstMetadata")
)

######################## Jaccard Evaluation ####################################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "beta_cluster_eval_jaccard_df",
    {
        perform_cluster_eval("jaccard", full_comp_df)
    },
    depends = c("mgtstMetadata")
)

######################## Unifrac Metric Evaluation #############################
## Execute evaluation and cache results 

ProjectTemplate::cache(
    "beta_cluster_eval_unifrac_df", 
    {
       perform_cluster_eval("unifrac", full_comp_df)
    }, 
    depends = c("mgtstMetadata")
)
