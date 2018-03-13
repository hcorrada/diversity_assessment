## Generate input data for abundance curve and prevalence 
calc_sparsity <- function(ps){ 
    
    ## Excluding no template controls
    non_ntc_samples <- sample_data(ps)$biosample_id != "NTC"
    ps <- prune_samples(non_ntc_samples, ps)
    
    ## generate count matrix
    mat <- as.matrix(otu_table(ps)) 
    
    ## calculate sparsity
    nentry <- length(mat)
    nzero <- sum(mat == 0)
    
    ## sparsity 
    nzero/nentry
}

## Identify number of singletons in each dataset
calc_singletons <- function(ps){ 
    
    ## Excluding no template controls
    non_ntc_samples <- sample_data(ps)$biosample_id != "NTC"
    ps <- prune_samples(non_ntc_samples, ps)
    
    ## calculate number of singletons
    nsingle <- length(which(taxa_sums(ps)==1))
    
    ## singletons
    nsingle
}

## Pipeline characterization
## Columns 
# - `pipe` - pipeline name  
# - `n_taxa` - number of features in count table  
# - `n_samples` - number of samples in count table  
# - `per_sample_info` - list column with data frames with per sample info
#   - `sample_id` - unique sample identifier, consistent `sample_id` in `mgtstMetadata` 
#   - `total_abu` - count table total per sample abundance  
#   - `Observed` - per sample number of observed features  
#   - `read` - per sample number of sequencing reads 
#   - `pass_rate` - proportion of reads passing quality filter  
 
get_pipe_info <- function(){
    
    pipe_info_df <- tibble(rds_file = list.files("data/phyloseq_objects", full.names = TRUE)) %>%
        mutate(
            pipe = str_extract(rds_file, "(?<=ts/).*(?=_ps)"),
            ps = map(rds_file, readRDS),
            n_taxa = map_dbl(ps, ntaxa),
            num_singletons = map_dbl(ps, calc_singletons),
            n_samples = map_dbl(ps, nsamples), 
            sparsity = map_dbl(ps, calc_sparsity), 
            taxa_per_sample = map(ps, estimate_richness, measures = "Observed"),
            taxa_per_sample = map(taxa_per_sample, rownames_to_column, var = "sample_id"),
            total_abu = map(ps, sample_sums)
        ) %>%
        select(-ps, -rds_file) 

    lib_size_df <- seq_char_df %>% 
        filter(read_dir == "R1") %>% 
        select(sample_id, read) %>% as_tibble()
    
    pipe_info_df %>% 
        unnest() %>% 
        mutate(sample_id = str_replace(sample_id, "\\.","-")) %>% 
        left_join(lib_size_df) %>% 
        mutate(pass_rate = total_abu / read) %>% 
        select(pipe, n_taxa, n_samples, sparsity, num_singletons, sample_id, 
               Observed, total_abu, read, pass_rate) %>% 
        group_by(pipe, n_taxa, n_samples, sparsity, num_singletons) %>% 
        nest() %>% 
        dplyr::rename(per_sample_info = data)
        
}

ProjectTemplate::cache(variable = "pipe_char_df", 
                       CODE = {get_pipe_info()}, 
                       depends = c("mgtstMetadata","seq_char_df"))


## Identify number of reads per sample at each rarefaction level
calc_rare_lev <- function(ps){ 
    
    ## calculate number of reads at rarefaction level
    rare_lev <- unique(sample_sums(ps))
    
    rare_lev
}

## Rarefaction characterization
## Columns 
# - `pipe` - pipeline name 
# - `normalization` - normalization method
# - `n_taxa` - number of features in count table  
# - `n_samples` - number of samples in count table  

get_rare_info <- function(){
    rds_file = list.files("data/norm_data", full.names = TRUE)
    rare_info_df <- tibble(file=rds_file[grep("rare", rds_file)]) %>%
        mutate(
            pipe = gsub(file, pattern="(data/norm_data/)([a-z].*)_(rare[0-9q][0-9]*).rds", replacement="\\2"),
            norm = gsub(file, pattern="(data/norm_data/)([a-z].*)_(rare[0-9q][0-9]*).rds", replacement="\\3"),
            ps = map(file, readRDS),
            n_reads = map_dbl(ps, calc_rare_lev),
            n_taxa = map_dbl(ps, ntaxa),
            n_samples = map_dbl(ps, nsamples)
           ) %>%
        select(-ps, -file) 
}

ProjectTemplate::cache(variable = "rare_char_df", 
                       CODE = {get_rare_info()})
