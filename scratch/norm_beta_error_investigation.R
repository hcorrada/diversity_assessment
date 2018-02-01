### Error Eval

## Normalization Results
norm_files <- list.files("data/norm_data", full.names = T) %>%
    set_names(basename(.))


norm_error <- norm_files %>% 
    map(readRDS) %>% 
    map(pluck, "error")

## No errors for pipeline normalization 
norm_error %>% keep(~!is.null(.))


## Beta diversity estimate results  
beta_div_files <- list.files("data/diversity_data", full.names = T) %>%
    set_names(basename(.)) 

beta_div_error <- beta_div_files %>% 
    map(readRDS) %>% 
    map(pluck, "error") 


## Error output includes error message and call, call contains matrix metric was calculated on resulting in large list, plucking message for error analysis
beta_div_error_info <- beta_div_error %>% 
    keep(~!is.null(.)) %>% 
    map(pluck, "message")

## Large list 1.9 Gb removing from memory
rm(beta_div_error)

## Unable to calculate diversity for 8 normalized datasets 
length(beta_div_error_info)

## All Bray Curtis estimates but different pipelines and normalization methods  
names(beta_div_error_info)

## [1] "dada_CSS_bray.rds"          "dada_TSS_bray.rds"          "dada_UQ_bray.rds"           "deblur_UQ_bray.rds" 
## [5] "mothur_UQ_bray.rds"         "qiimeClosedRef_UQ_bray.rds" "qiimeDeNovo_UQ_bray.rds"    "qiimeOpenRef_UQ_bray.rds"  

## -- no errors after fixing normalization bugs

# Error output contains call as well as error message only looking at error message
beta_div_error_info 

## All NA/NaN/Inf in foreign function call 
dada_CSS <- readRDS("data/norm_data/dada_CSS.rds")
count_dat <- otu_table(dada_CSS)

## Featues with All NA
bad_features <- colnames(count_dat)[is.na(colSums(count_dat))]
length(bad_features)

## Check ps object for NAs in count table   
dada_ps <- readRDS("data/phyloseq_objects/dada_ps.rds")
bad_features_ps <- phyloseq::prune_taxa(bad_features, dada_ps)
colSums(otu_table(bad_features_ps))
non_zero_samples <- prune_samples(rowSums(otu_table(bad_features_ps)) > 0, bad_features_ps)
bad_nzero_otu_table <- otu_table(non_zero_samples)

## Removing no template controls
remove_ntc <- function(ps){
non_ntc_samples <- sample_data(ps)$biosample_id != "NTC"
prune_samples(non_ntc_samples, ps)
}
## Removing samples with no reads
remove_no_read_samples <- function(ps) prune_samples(sample_sums(ps) > 0,  ps)
dada_no_ntc <- dada_ps %>% remove_ntc() %>% remove_no_read_samples()
ps <- dada_no_ntc
method <- "CSS"
require(matrixStats)
## Extract count matrix from phyloseq object
count_mat <- as(otu_table(ps), "matrix")
## Cumulative sum scaling Paulson et al. 2013
norm_factors <- metagenomeSeq::calcNormFactors(as(count_mat, "matrix"), p = 0.75)
norm_factors <- norm_factors$normFactors
## Normalizing counts
norm_mat <- sweep(count_mat, 2, norm_factors,'/')
norm_mat
dims(norm_mat)
dim(norm_mat)
colnames(norm_mat)
norm_mat[colnames(norm_mat) %in% bad_features,]
norm_mat[,colnames(norm_mat) %in% bad_features]

## NA norm factors 
norm_factors
taxa_are_rows(ps)
taxa_are_rows(dada_ps)
## Extract count matrix from phyloseq object
count_mat <- as(otu_table(ps), "matrix")
## Consistent format - taxa as rows
if (!taxa_are_rows(ps)) {
count_mat <- t(count_mat)
}
## Cumulative sum scaling Paulson et al. 2013
norm_factors <- metagenomeSeq::calcNormFactors(count_mat, p = 0.75)
norm_factors <- norm_factors$normFactors
norm_factors
## Normalizing counts
norm_mat <- sweep(count_mat, 2, norm_factors,'/')
norm_mat[rownames(norm_mat) %in% bad_features,]

## Need to rerun normalization and beta diversity phyloseq objects where taxa are cols  
get_taxa_are_rows <- compose(readRDS, taxa_are_rows)

## Only DADA2 has taxa as rows
list.files("data/phyloseq_objects", full.names = TRUE) %>% 
    map(readRDS) %>% 
    map(taxa_are_rows)

## Need to figure out issues with other failed runs 


## Mothur UQ  
## All NA/NaN/Inf in foreign function call 
mothur_UQ <- readRDS("data/norm_data/mothur_UQ.rds")
count_dat <- otu_table(mothur_UQ)

## Featues with All NA
bad_features <- rownames(count_dat)[is.na(rowSums(count_dat))]
length(bad_features)

## Issue with normalization function - UQ sets 0's to NA, needed to reset to 0 after getting norm factors
prune_taxa(bad_features[1:10], mothur_UQ) %>% otu_table()
