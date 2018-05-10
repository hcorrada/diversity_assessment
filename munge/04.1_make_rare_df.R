## Data frame for rarefaction curve plot

ps_iNEXT <- function(ps){
    # sample_set <- c("nist_run1_2-A10","nist_run2_2-A10",
    #                 "jhu_run1_2-A10","jhu_run2_2-A10")
    seq_lab <- sample_data(ps)$seq_lab 
    seq_run <- sample_data(ps)$seq_run
    sample_data(ps)$seq_lab_run <- paste0(seq_lab, seq_run)
    ps <- merge_samples(ps, group = "seq_lab_run")
    
    
    
    count_tbl <- ps %>% 
        # prune_samples(sample_set, ps)  %>%
    {prune_taxa(taxa_sums(.) > 0, .)} %>%
        otu_table()
    
    if (!taxa_are_rows(ps)) {
        count_tbl <- t(count_tbl)
    }
    
    count_df <- as.data.frame(count_tbl)
    
    iNEXT::iNEXT(count_df)
}


make_rare_plot_df <- function(){
    ps_list <- list.files("data/phyloseq_objects", full.names = TRUE) %>%
        set_names(str_remove(basename(.), "_ps.rds")) %>%
        map(readRDS)
    
    tibble(ps_obj = ps_list) %>%
        add_column(pipe = names(ps_list)) %>%
        mutate(pipe = case_when(pipe == "deblur" ~ "q_deblur",
                                pipe == "qiimeOpenRef" ~ "q_open",
                                pipe == "qiimeClosedRef" ~ "q_closed",
                                pipe == "qiimeDeNovo"~ "q_denovo",
                                TRUE ~ pipe),
               pipe = factor(pipe)) %>%
        mutate(inext_dat = map(ps_obj, ps_iNEXT),
               rare_plot = map2(inext_dat, pipe, ~{iNEXT::ggiNEXT(.x) +
                       labs(x = "Sampling Depth", y = "Feature Diversity") +
                       theme_bw() +
                       ggtitle(.y) +
                       scale_shape_manual(values = run_shapes) +
                       scale_color_manual(values=run_colors2)})) %>%
        dplyr::select(-ps_obj, -inext_dat)
}



## Generate and cache data frame 
ProjectTemplate::cache("rare_plot_df", 
                       {make_rare_plot_df()},
                       depends = c("mgtstMetadata"))