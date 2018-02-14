# setwd("~/diversity_assessment/")
library(ProjectTemplate)
load.project()


factors <- c("mean_observed_feat", "cov_observed_feat", 
           "mean_total_abu", "cov_total_abu", 
           "mean_num_reads", "cov_num_reads", 
           "mean_pass_rate", "cov_pass_rate")

factor_labels <- c("Mean Observed Features", "Coefficient of Variation- Observed Features", 
           "Mean Total Abundance", "Coefficient of Variation- Total Abundance", 
           "Mean Number of Reads", "Coefficient of Variation- Number of Reads", 
           "Mean Pass Rate", "Coefficient of Variation- Pass Rate")

norm_methods <- c("RAW", "rare2000", "rare5000", "rare10000", "rareq15", 
                  "CSS", "TMM", "RLE", "TSS","UQ")

seq_char_comparisons$normalization <- factor(seq_char_comparisons$normalization, 
                                           levels = norm_methods,
                                           ordered=T)

pdf("reports/bdiv_sequencing_characteristics.pdf", height=8, width=12)
for(i in 1:length(factors)){
    
    for(j in c("jaccard", "unifrac", "wunifrac", "bray")){
        a<-ggplot(subset(seq_char_comparisons, metric==j), 
                  aes_string(x=factors[i], y="mean_dist"))+
            geom_point()+geom_smooth(method = "lm")+
            facet_grid(pipe~normalization, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(factor_labels[i])+ylab("Mean Distance")+ ggtitle(j)
        print(a)
        
    }
    
    b<-ggplot(subset(seq_char_comparisons, normalization %in% c("RAW", "rare2000", "rare5000", "rare10000", "rareq15")), 
              aes_string(x=factors[i], y="mean_dist", 
                         group="pipe", color="pipe", 
                         fill="pipe"))+
        geom_smooth(method = "lm")+facet_grid(metric~normalization, scales="free")+theme_bw()+
        theme(axis.text.x=element_text(angle = -45, hjust = 0))+
        xlab(factor_labels[i])+ylab("Mean Distance")
    
    c<-ggplot(subset(seq_char_comparisons, normalization %in% c("RAW", "rare2000", "rare5000", "rare10000", "rareq15")), 
           aes_string(x=factors[i], y="mean_dist", 
               group="normalization", color="normalization", 
               fill="normalization"))+
        geom_smooth(method = "lm")+facet_grid(metric~pipe, scales="free")+theme_bw()+
        theme(axis.text.x=element_text(angle = -45, hjust = 0))+
        xlab(factor_labels[i])+ylab("Mean Distance")
    
    d<-ggplot(subset(seq_char_comparisons, normalization %in% c("RAW", "CSS", "TMM", "RLE", "TSS","UQ") & metric %in% c("bray", "wunifrac")), 
              aes_string(x=factors[i], y="mean_dist", 
                         group="pipe", color="pipe", 
                         fill="pipe"))+
        geom_smooth(method = "lm")+facet_grid(metric~normalization, scales="free")+theme_bw()+
        theme(axis.text.x=element_text(angle = -45, hjust = 0))+
        xlab(factor_labels[i])+ylab("Mean Distance")
    
    e<-ggplot(subset(seq_char_comparisons, normalization %in% c("RAW", "CSS", "TMM", "RLE", "TSS","UQ") & metric %in% c("bray", "wunifrac")), 
              aes_string(x=factors[i], y="mean_dist", 
                         group="normalization", color="normalization", 
                         fill="normalization"))+
        geom_smooth(method = "lm")+facet_grid(metric~pipe, scales="free")+theme_bw()+
        theme(axis.text.x=element_text(angle = -45, hjust = 0))+
        xlab(factor_labels[i])+ylab("Mean Distance")
    
    print(b)
    print(c)
    print(d)
    print(e)
}
dev.off()

factors<-c("mean_quality", "cov_quality")

factor_labels<-c("Mean Quality", "Coefficient of Variation- Quality")

seq_qual_comparisons$normalization<-factor(seq_qual_comparisons$normalization, 
                                           levels= c("RAW", "rare2000", "rare5000", "rare10000", "rareq15", 
                                                     "CSS", "TMM", "RLE", "TSS","UQ"),
                                           ordered=T)


pdf("reports/bdiv_sequencing_quality.pdf", height=8, width=12)
for(i in 1:length(factors)){

    for(j in c("jaccard", "unifrac", "wunifrac", "bray")){
        
        a<-ggplot(subset(seq_qual_comparisons, metric==j), 
                  aes_string(x=factors[i], y="mean_dist", 
                             color="read_dir", fill="read_dir",
                             group="read_dir"))+
            geom_point()+geom_smooth(method = "lm")+
            facet_grid(pipe~normalization, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(factor_labels[i])+ylab("Mean Distance") + ggtitle(j)
        
        print(a)
        
    }
    
    for(direction in c("R1", "R2")){
        
        b<-ggplot(subset(seq_qual_comparisons, read_dir==direction & normalization %in% c("RAW", "rare2000", "rare5000", "rare10000", "rareq15")), 
                  aes_string(x=factors[i], y="mean_dist", 
                             group="pipe", color="pipe", 
                             fill="pipe"))+
            geom_smooth(method = "lm")+facet_grid(metric~normalization, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(paste0(factor_labels[i],"_", direction))+ylab("Mean Distance")
        
        c<-ggplot(subset(seq_qual_comparisons, read_dir==direction & normalization %in% c("RAW", "rare2000", "rare5000", "rare10000", "rareq15")), 
                  aes_string(x=factors[i], y="mean_dist", 
                             group="normalization", color="normalization", 
                             fill="normalization"))+
            geom_smooth(method = "lm")+facet_grid(metric~pipe, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(paste0(factor_labels[i],"_", direction))+ylab("Mean Distance")
        
        d<-ggplot(subset(seq_qual_comparisons, read_dir==direction & normalization %in% c("RAW", "CSS", "TMM", "RLE", "TSS","UQ") & metric %in% c("bray", "wunifrac")), 
                  aes_string(x=factors[i], y="mean_dist", 
                             group="pipe", color="pipe", 
                             fill="pipe"))+
            geom_smooth(method = "lm")+facet_grid(metric~normalization, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(paste0(factor_labels[i],"_", direction))+ylab("Mean Distance")
        
        e<-ggplot(subset(seq_qual_comparisons, read_dir==direction & normalization %in% c("RAW", "CSS", "TMM", "RLE", "TSS","UQ") & metric %in% c("bray", "wunifrac")), 
                  aes_string(x=factors[i], y="mean_dist", 
                             group="normalization", color="normalization", 
                             fill="normalization"))+
            geom_smooth(method = "lm")+facet_grid(metric~pipe, scales="free")+theme_bw()+
            theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            xlab(paste0(factor_labels[i],"_", direction))+ylab("Mean Distance")
        
        print(b)
        print(c)
        print(d)
        print(e)
    }
    
}
dev.off()
