## Evaluation of biological v. technical replicates
## beta diversity metrics used in McMurdie and Holmes 2014 and Weiss et al 2017
## comparing unmixed pre-exposure samples to titration and post-exposure samples

##################### Functions for PAM Cluster Evaluation #####################

#' Evaluate biological versus technical variation
#'
#' @param map mapping file with sample information 
#' @param metric diversity metric of interest (bray, jaccard, unifrac or wunifrac)
#'
#' @return data frame with diversity value for paired samples for all pipelines and normalization methods
#' @export
#'
#' @examples
compute_diversity_comparisons_with_seq_char<-function(metric, sample_comparisons){
    
    # Generate list of diversity matrices for piplines/norm methods
    data<-make_beta_div_df(metric)
    
    # For each beta div matrix
    beta_div_comp<- lapply(1:length(data$pipe), function(i){
        # Extract matrix
        beta_div<-as.data.frame(as.matrix(data$dist_results[[i]][[1]]))
        # Keep only matrix with samples of interest
        beta_div_m<-as.data.frame(beta_div)[which(row.names(beta_div) %in% as.character(unique(c(sample_comparisons$sample_a, sample_comparisons$sample_b)))), 
                                            which(colnames(beta_div) %in% as.character(unique(c(sample_comparisons$sample_a, sample_comparisons$sample_b))))]
        # Melt dataframe so we can look at matched samples 
        beta_div_m$sample<-row.names(beta_div_m)
        beta_div_m2<-gather(beta_div_m, sample)
        colnames(beta_div_m2)<-c("sample_a", "sample_b", "value")
        # Look only at matched samples
        beta_div_m2<-as.data.frame(merge(beta_div_m2, sample_comparisons))
        beta_div_m2$value<-as.numeric(as.character(beta_div_m2$value))
        # Calculate average distance for pcr replicates
        beta_div_summary<- beta_div_m2 %>% 
            group_by(replicate) %>%
            summarise(mean_dist=mean(value), N=length(value))
        
        # Calculate sequencing info for replicates
        map<-mgtstMetadata
        map$replicate<-paste0(map$biosample_id, "_t", 
                              map$t_fctr, "_",
                              map$seq_lab, "_run",
                              map$seq_run)
        map<-map[,c("sample_id", "replicate")]
        seq_info<-pipe_char_df$per_sample_info[[which(pipe_char_df$pipe == data$pipe[i])]]
        seq_info<-merge(map, seq_info)
        seq_info_summary<- seq_info %>% 
            group_by(replicate) %>%
            summarise(mean_observed_feat=mean(Observed),
                      sd_observed_feat=sd(Observed),
                      cov_observed_feat=sd_observed_feat/mean_observed_feat,
                      mean_total_abu=mean(total_abu),
                      sd_total_abu=sd(total_abu),
                      cov_total_abu=sd_total_abu/mean_total_abu,
                      mean_num_reads=mean(read),
                      sd_num_reads=sd(read),
                      cov_num_reads=sd_num_reads/mean_num_reads,
                      mean_pass_rate=mean(pass_rate),
                      sd_pass_rate=sd(pass_rate),
                      cov_pass_rate=sd_pass_rate/mean_pass_rate,
                      N_seq_info=length(sample_id))
        
        # Add additional info to melted dataframe
        summary<-merge(beta_div_summary, seq_info_summary)
        summary$pipe<-data$pipe[i]
        summary$normalization<-data$method[i]
        summary$metric<-metric
        
        return(summary)
    })
    # Merge diversity tables from all pipelines
    output<-do.call("rbind", beta_div_comp)
    return(output)
    #return(beta_div_comp)
}


######## Extract biol or tech variation from diversity matrices ######################

# Identify comparisons within PCR replicates
# Look at all combinations of samples
sample_comparisons<-as.data.frame(t(combn(as.character(mgtstMetadata$sample_id), 2)))
colnames(sample_comparisons)<-c("sample_a", "sample_b")

# Determine what characteristics are shared by each sample pair
sample_comparisons<-merge(sample_comparisons, mgtstMetadata, by.x="sample_a", by.y="sample_id")
colnames(sample_comparisons)[3:ncol(sample_comparisons)]<-
    paste0( "sample_a_",colnames(sample_comparisons)[3:ncol(sample_comparisons)])
sample_comparisons<-merge(sample_comparisons, mgtstMetadata, by.x="sample_b", by.y="sample_id")
colnames(sample_comparisons)[9:ncol(sample_comparisons)]<-
    paste0("sample_b_",colnames(sample_comparisons)[9:ncol(sample_comparisons)])

sample_comparisons<-subset(sample_comparisons, 
                           sample_a_biosample_id == sample_b_biosample_id &
                               sample_a_t_fctr == sample_b_t_fctr &
                               sample_a_seq_lab == sample_b_seq_lab & 
                               sample_a_seq_run == sample_b_seq_run)


sample_comparisons$replicate<-paste0(sample_comparisons$sample_a_biosample_id, "_t", 
                                     sample_comparisons$sample_a_t_fctr, "_",
                                     sample_comparisons$sample_a_seq_lab, "_run",
                                     sample_comparisons$sample_a_seq_run)

sample_comparisons<-sample_comparisons[,c("sample_a", "sample_b", "replicate")]

# Unweighted metrics
unweighted_unifrac_comp<-compute_diversity_comparisons_with_seq_char("unifrac", sample_comparisons)

jaccard_comp<-compute_diversity_comparisons_with_seq_char("jaccard", sample_comparisons)

# Weighted metrics
weighted_unifrac_comp<-compute_diversity_comparisons_with_seq_char("wunifrac", sample_comparisons)

bray_curtis_comp<-compute_diversity_comparisons_with_seq_char("bray", sample_comparisons)

seq_char_comparisons<-rbind(unweighted_unifrac_comp, jaccard_comp)
seq_char_comparisons<-rbind(seq_char_comparisons, weighted_unifrac_comp)
seq_char_comparisons<-rbind(seq_char_comparisons, bray_curtis_comp)

# Only keep comparisons where you have at least three of four replicates
seq_char_comparisons<-subset(seq_char_comparisons, N_seq_info >=3 & N >=3)
# IF ONLY 3 SAMPLES-- SOMETIMES LOOKING AT QUALITY OF FOUR SAMPLES COMBINED-- 
# NEED TO GO BACK AND FIX!!!
ProjectTemplate::cache('seq_char_comparisons', 
                       depends = c("mgtstMetadata"))