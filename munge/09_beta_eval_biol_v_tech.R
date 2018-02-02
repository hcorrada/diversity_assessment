## Evaluation of biological v. technical replicates
## beta diversity metrics used in McMurdie and Holmes 2014 and Weiss et al 2017
## comparing unmixed pre-exposure samples to titration and post-exposure samples

##################### Functions for Biol. v. Tech. Replicate Evaluation #####################

#' Evaluate biological versus technical variation
#'
#' @param map mapping file with sample information 
#' @param metric diversity metric of interest (bray, jaccard, unifrac or wunifrac)
#' @param variation_tests dataframe containing which samples should be compared for different parameters
#' 
#' @return data frame with diversity value for paired samples for all pipelines and normalization methods
#' @export
#'
#' @examples
compute_diversity_comparisons<-function(map, metric, variation_tests){
    
    # Generate list of diversity matrices for piplines/norm methods
    data<-make_beta_div_df(metric)
    
    # For each beta div matrix
    beta_div_comp<- lapply(1:length(data$pipe), function(i){
        # Extract matrix
        beta_div<-as.matrix(data$dist_results[[i]][[1]])
        # Keep only matrix with samples of interest
        beta_div_m<-as.data.frame(beta_div)[which(row.names(beta_div) %in% as.character(map$sample_id)), 
                                            which(colnames(beta_div) %in% as.character(map$sample_id))]
        # Melt dataframe so we can look at matched samples 
        beta_div_m$sample<-row.names(beta_div_m)
        beta_div_m2<-gather(beta_div_m, sample)
        colnames(beta_div_m2)<-c("sample_a", "sample_b", "value")
        # Add additional info to melted dataframe
        beta_div_m2$pipe<-data$pipe[i]
        beta_div_m2$normalization<-data$method[i]
        beta_div_m2$metric<-metric
        return(beta_div_m2)
    })
    # Merge diversity tables from all pipelines
    output<-do.call("rbind", beta_div_comp)
    # Merge with predetermined diversity tests
    output2<-merge(variation_tests, output, by=c("sample_a", "sample_b"))
    output2$normalization<-factor(output2$normalization, levels=c("RAW",
                                                                  "rare2000",
                                                                  "rare5000",
                                                                  "rare10000",
                                                                  "rareq15",
                                                                  "CSS",
                                                                  "RLE",
                                                                  "TMM",
                                                                  "TSS",
                                                                  "UQ"),
                                  ordered = T)
    output2$normalization_type<-NA
    output2$normalization_type[which(output2$normalization == "RAW")]<-c("none")
    output2$normalization_type[which(output2$normalization %in% c("rare2000", 
                                                                  "rare5000",
                                                                  "rare10000",
                                                                  "rareq15"))]<-c("rarefaction")
    output2$normalization_type[which(output2$normalization %in% c("CSS", 
                                                                  "RLE",
                                                                  "TMM",
                                                                  "TSS",
                                                                  "UQ"))]<-c("abundance_based")
    return(data.frame(output2))
}


######## Identify paired samples demonstrating biol or tech variation ######################

# Extract on only the unmixed pre samples (t_fctr==20) or unmixed post samples (t_fctr==0)
map_sub<-subset(mgtstMetadata, 
                mgtstMetadata$t_fctr %in% c(0, 20) & mgtstMetadata$biosample_id != "NTC")
# Fix factoring of some columns
map_sub<-droplevels(map_sub)

# Look at all combinations of samples
sample_comparisons<-as.data.frame(t(combn(as.character(map_sub$sample_id), 2)))
colnames(sample_comparisons)<-c("sample_a", "sample_b")

# Determine what characteristics are shared by each sample pair
sample_comparisons<-merge(sample_comparisons, map_sub, by.x="sample_a", by.y="sample_id")
colnames(sample_comparisons)[3:ncol(sample_comparisons)]<-
    paste0( "sample_a_",colnames(sample_comparisons)[3:ncol(sample_comparisons)])
sample_comparisons<-merge(sample_comparisons, map_sub, by.x="sample_b", by.y="sample_id")
colnames(sample_comparisons)[9:ncol(sample_comparisons)]<-
    paste0("sample_b_",colnames(sample_comparisons)[9:ncol(sample_comparisons)])

# Same individual
sample_comparisons$same_individual<-
    (sample_comparisons$sample_a_biosample_id==sample_comparisons$sample_b_biosample_id)

# Same timepoint
sample_comparisons$same_timepoint<-
    (sample_comparisons$sample_a_t_fctr==sample_comparisons$sample_b_t_fctr)

# Same PCR product 
sample_comparisons$same_pcr_product<-
    ((sample_comparisons$sample_a_pcr_16S_plate==sample_comparisons$sample_b_pcr_16S_plate) & 
         (sample_comparisons$sample_a_pos == sample_comparisons$sample_b_pos))

# Same sequencing lab
sample_comparisons$same_seq_lab<-
    (sample_comparisons$sample_a_seq_lab==sample_comparisons$sample_b_seq_lab)

# Same sequencing run
sample_comparisons$same_seq_run<-
    (sample_comparisons$sample_a_seq_run==sample_comparisons$sample_b_seq_run)

# Remove unused columns
sample_comparisons<-sample_comparisons[,-which(colnames(sample_comparisons) %in% 
                                                   c("sample_a_biosample_id", "sample_b_biosample_id",
                                                     "sample_a_t_fctr", "sample_b_t_fctr",
                                                     "sample_a_pcr_16S_plate", "sample_b_pcr_16S_plate",
                                                     "sample_a_pos", "sample_b_pos",
                                                     "sample_a_seq_lab", "sample_b_seq_lab",
                                                     "sample_a_seq_run", "sample_b_seq_run"))]

# Classify variation (biological or technical) and add variation labels (between timepoint variation, etc.)
## biological variation
# between timepoint variation
between_timepoint<-subset(sample_comparisons, sample_comparisons$same_timepoint==FALSE &
                              sample_comparisons$same_seq_lab==TRUE &
                              sample_comparisons$same_seq_run==TRUE)
between_timepoint$variation<-"biological"
between_timepoint$variation_label<-"btw_time"

# between individual variation
between_individual<-subset(sample_comparisons, sample_comparisons$same_individual==FALSE & 
                               sample_comparisons$same_timepoint==TRUE &
                               sample_comparisons$same_seq_lab==TRUE &
                               sample_comparisons$same_seq_run==TRUE)
between_individual$variation<-"biological"
between_individual$variation_label<-"btw_subj_w/in_time"

# within individual timepoint variation
within_individual_timepoint<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                                        sample_comparisons$same_timepoint==FALSE &
                                        sample_comparisons$same_seq_lab==TRUE &
                                        sample_comparisons$same_seq_run==TRUE)
within_individual_timepoint$variation<-"biological"
within_individual_timepoint$variation_label<-"w/in_subj_btw_time"

variation<-rbind(between_timepoint, between_individual)
variation<-rbind(variation, within_individual_timepoint)
rm(between_timepoint, between_individual, within_individual_timepoint)

## technical variation
# between seq labs
between_seq_labs<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                             sample_comparisons$same_timepoint==TRUE &
                             sample_comparisons$same_pcr_product==TRUE &
                             sample_comparisons$same_seq_lab==FALSE)
between_seq_labs$variation<-"technical"
between_seq_labs$variation_label<-"btw_labs"

# within seq labs run variation
within_seq_lab_runs<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                                sample_comparisons$same_timepoint==TRUE &
                                sample_comparisons$same_seq_lab==TRUE &
                                sample_comparisons$same_pcr_product==TRUE &
                                sample_comparisons$same_seq_run==FALSE)
within_seq_lab_runs$variation<-"technical"
within_seq_lab_runs$variation_label<-"w/in_lab_runs"

# within seq labs PCR product variation
within_seq_lab_pcr<-subset(sample_comparisons, sample_comparisons$same_individual==TRUE & 
                               sample_comparisons$same_timepoint==TRUE &
                               sample_comparisons$same_seq_lab==TRUE &
                               sample_comparisons$same_seq_run==TRUE &
                               sample_comparisons$same_pcr_product==FALSE)
within_seq_lab_pcr$variation<-"technical"
within_seq_lab_pcr$variation_label<-"w/in_lab_pcr"

variation<-rbind(variation, between_seq_labs)
variation<-rbind(variation, within_seq_lab_runs)
variation<-rbind(variation, within_seq_lab_pcr)
rm(between_seq_labs, within_seq_lab_runs, within_seq_lab_pcr)

# Order variation types
variation$variation_label<-factor(variation$variation_label, 
                                  levels=c("btw_subj_w/in_time", 
                                           "w/in_subj_btw_time", 
                                           "btw_time", "w/in_lab_pcr", 
                                           "w/in_lab_runs","btw_labs"))

# Generate dataframe to plot what comparisons are being made
variation_tmp<-unique(variation[,c(3:7,9)])
# Clean up and merge duplicates in same variation label
variation_tmp[1,1]<-"TRUE OR FALSE"
variation_tmp<-variation_tmp[-2,]

variation_tmp[4,5]<-"TRUE OR FALSE"
variation_tmp<-variation_tmp[-5,]

biol_v_tech_variation_comparison_map<-gather(variation_tmp, variation_label, variable)
rm(variation_tmp)

biol_v_tech_variation_comparison_map$value<-factor(biol_v_tech_variation_comparison_map$value, levels = c("TRUE", "TRUE OR FALSE", "FALSE"), ordered = TRUE)

######## Extract biol or tech variation from diversity matrices ######################

# Unweighted metrics
unweighted_unifrac_comp<-compute_diversity_comparisons(map_sub, "unifrac", variation)
jaccard_comp<-compute_diversity_comparisons(map_sub, "jaccard", variation)

# Weighted metrics
weighted_unifrac_comp<-compute_diversity_comparisons(map_sub, "wunifrac", variation)
bray_curtis_comp<-compute_diversity_comparisons(map_sub, "bray", variation)

# Only need to calculate for one metric since it should be same number of overall comparisons for all
biol_v_tech_variation_comparison_numbers<- weighted_unifrac_comp %>% 
    group_by(variation, variation_label, pipe, normalization, normalization_type) %>% 
    summarise(number_of_comparisons=length(variation_label)) 

biol_v_tech_variation_comparisons<-rbind(bray_curtis_comp, jaccard_comp)
biol_v_tech_variation_comparisons<-rbind(biol_v_tech_variation_comparisons, unweighted_unifrac_comp)
biol_v_tech_variation_comparisons<-rbind(biol_v_tech_variation_comparisons, weighted_unifrac_comp)

biol_v_tech_variation=list(biol_v_tech_variation_comparison_map, 
                           biol_v_tech_variation_comparison_numbers, 
                           biol_v_tech_variation_comparisons)
########################  Cache results ################################

ProjectTemplate::cache('biol_v_tech_variation', 
                       depends = c("mgtstMetadata"))

