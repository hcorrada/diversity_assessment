# ## Generating Data Frame with weighted and unweighted diversity metrics
# diversity_rds <-
#     list.files("data/diversity_data", full.names = TRUE)
# 
# diversity_names <-
#     basename(diversity_rds) %>% str_replace(".rds", "")
# 
# diversity_df <- data_frame(div_info = diversity_names,
#                            rds_file = diversity_rds) %>%
#     separate(div_info, c("pipe", "method", "dist_method")) %>%
#     mutate(dist_method = paste0(dist_method, "_dist"))
# 
# ############## Unweighted metrics - rarified data only #########################
# ## Unweighted metrics include unweighted unifrac and jaccard
# ProjectTemplate::cache("unweighted_beta_df", {
#     diversity_df %>%
#         filter(dist_method %in% c("unifrac_dist", "jaccard_dist")) %>%
#         mutate(dist_results = map(rds_file, readRDS)) %>%
#         select(-rds_file)
# })
# 
# ########### Weighted Metrics ###################################################
# ## Weighted metrics evaluated compared include weighted unifrac and bray-curtis
# ProjectTemplate::cache("weighted_beta_df", {
#     diversity_df %>%
#         filter(dist_method %in% c("wunifrac_dist", "bray_dist")) %>%
#         mutate(dist_results = map(rds_file, readRDS)) %>%
#         select(-rds_file)
# })