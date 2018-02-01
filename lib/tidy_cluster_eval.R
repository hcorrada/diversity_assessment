## Tidy cluster evaluation results - not returning tidy objects for debugging
tidy_cluster_eval <- function(cluster_eval_df) {
  cluster_eval_df %>%
      mutate(eval_result = map(eval_output, pluck, "result")) %>%
      mutate(eval_error = map_lgl(eval_result, is.null)) %>%
      filter(!eval_error) %>%
      select(pipe, method, dist_method, eval_result) %>%
      unnest() %>%
      select(-comp_df, -cluster_assignments)
}