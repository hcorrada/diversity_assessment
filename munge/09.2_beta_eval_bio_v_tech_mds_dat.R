## Data for PCoA ordination sanity check plot

## Function to prepare phyloseq objects for analysis 

get_ord_df <- function(){ 
    ## Normalized phyloseq objects
    dada_tss_full <- readRDS("data/norm_data/dada_TSS.rds")
    dada_tmm_full <- readRDS("data/norm_data/dada_TMM.rds")
    mothur_tss_full <- readRDS("data/norm_data/mothur_TSS.rds")
    mothur_tmm_full <- readRDS("data/norm_data/mothur_TMM.rds")
    
    ## Helper functions 
    subset_ps <- function(ps){
        ps <- subset_samples(ps, t_fctr %in% c(0)) 
        filter_taxa(ps, function(x) sum(x) > 0, TRUE)
    }

    factor_var <- function(ps){
        sample_data(ps)$t_fctr <- factor(sample_data(ps)$t_fctr)
        sample_data(ps)$seq_run <- paste0(sample_data(ps)$seq_lab, 
                                      sample_data(ps)$seq_run)
    
        ps
    }
    
    prep_ps <- purrr::compose(subset_ps, factor_var)
    
    dada_ps_df <- tibble(pipe = "DADA2", 
                         normalization = c("TSS", "TMM")) %>% 
        add_column(ps = list(dada_tss_full, dada_tmm_full))
    
    mothur_ps_df <- tibble(pipe = "Mothur", 
                           normalization = c("TSS", "TMM")) %>% 
        add_column(ps = list(mothur_tss_full, mothur_tmm_full))
    
    make_ord_plot <- function(ps, ord, run_colors){
        plot_ordination(ps, ord, color = "seq_run", shape = "biosample_id") + 
            theme_bw()
    }
    
    ord_df <- bind_rows(dada_ps_df, mothur_ps_df) %>% 
        mutate(ps = map(ps, prep_ps)) %>% 
        mutate(Bray = map(ps, ordinate, method = "MDS"),
               `Weighted UniFrac` = map(ps, ordinate, method = "MDS", 
                                        distance = "wunifrac")) %>% 
        gather("dist_method", "ord", Bray, `Weighted UniFrac`)
    
    ord_plot_df <- ord_df %>%     
        mutate(ord_plot = map2(ps, ord, make_ord_plot)) %>% 
        arrange(pipe, normalization, dist_method)  
    
    ## Preparing data for plot
    ord_plot_df %>% 
        select(-ps, -ord) %>% 
        mutate(plot_dat = map(ord_plot, ~.$data),
               y_label = map_chr(ord_plot, ~.$labels$y),
               x_label = map_chr(ord_plot, ~.$labels$x)) %>% 
        mutate(x_label = str_extract(x_label, "(?<=\\[).*(?=\\])")) %>% 
        mutate(y_label = str_extract(y_label, "(?<=\\[).*(?=\\])")) %>% 
        mutate(pcoa_percent = paste0("Axis.1: ", x_label, "\nAxis.2: ", y_label)) %>% 
        select(-ord_plot) %>% 
        unnest()
}


## Generate and cache data frame 
ProjectTemplate::cache("ord_plot_dat_df", 
                       {get_ord_df()},
                       depends = c("mgtstMetadata"))
