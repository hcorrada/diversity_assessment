# diversity_assessment

## Data Descriptions  
* `seq_char_df` - tibble with sequence data characteristics  
    - `sample_id` - unique id for each sequenced sample   
    - `reads`` - number of reads  
    - `read_dir` - read direction   
    - `quality` - mode quality score  
    - `density` - proportion of reads with mode quality score 

* `pipe_char_df` - tibble with pipeline characteristics
    - `pipe` - pipeline name  
    - `n_taxa` - number of features in count table  
    - `n_samples` - number of samples in count table  
    - `per_sample_info` - list column with data frames with per sample info  
        - `sample_id` - unique sample identifier, consistent `sample_id` in `mgtstMetadata`  
        - `total_abu` - count table total per sample abundance   
        - `Observed` - per sample number of observed features    
        - `read` - per sample number of sequencing reads  
        - `pass_rate` - proportion of reads passing quality filter  

No longer in memory - too large ...
* `weighted_beta_df` and `unweighted_beta_df` are data frames with beta diversity estimates for the different pipelines calculated for each of the normalization methods.   
* The beta diversity distance metrics are included in the data frames as a list column with the distance matrix. 

## 
