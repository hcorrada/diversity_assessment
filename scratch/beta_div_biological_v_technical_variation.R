setwd("~/diversity_assessment/")
library(ProjectTemplate)
load.project()


# Plotting comparisons generated in munge/09_beta_eval_biol_v_tech and 
# saved in list biol_v_tech_variation

pdf("./reports/biological_v_technical_variation_diversity_plots.pdf", height=8, width=12)

## biol_v_tech_variation[[1]] is a dataframe containing information what features
## define bioligical replicates and technical replicates
print(ggplot(biol_v_tech_variation[[1]], aes(y=variation_label, x=variable, fill=value))+
  geom_tile()+theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  scale_fill_manual(values = c("black", "gray", "white"))+xlab("")+ylab("Type of Variation"))

## biol_v_tech_variation[[2]] is a dataframe containing information on how many comparsions
## are made per type of variation, pipeline, and normalization method tested

## Note: only weighted metrics were normalized using abundance based metrics

print(ggplot(subset(biol_v_tech_variation[[2]], normalization_type %in% c("none", "rarefaction")),
       aes(x=normalization, y=number_of_comparisons, fill=variation))+
  geom_bar(stat="identity")+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
  facet_grid(variation_label~pipe, scales="free_y")+xlab("Type of Variation")+
  scale_fill_manual(values = c("#d53e4f", "#3288bd")))

print(ggplot(subset(biol_v_tech_variation[[2]], normalization_type %in% c("none", "abundance_based")),
       aes(x=normalization, y=number_of_comparisons, fill=variation))+
    geom_bar(stat="identity")+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
    facet_grid(variation_label~pipe, scales="free_y")+xlab("Type of Variation")+
    scale_fill_manual(values = c("#d53e4f", "#3288bd")))

## biol_v_tech_variation[[3]] is a dataframe containing diversity values 
## for pairs of samples, per variation type, pipeline, and norm. methods
## for all diversity metrics tested

for(div_metric in unique(biol_v_tech_variation[[3]]$metric)){
    
    if(div_metric %in% c("wunifrac", "bray")){
        
        p<-ggplot(subset(biol_v_tech_variation[[3]], metric==div_metric & 
                         normalization_type %in% c("none","rarefaction")), 
                  aes(x=variation_label, y=value, fill=variation))+
            geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(div_metric)+
            scale_fill_manual(values = c("#d53e4f", "#3288bd"))
        
        q<-ggplot(subset(biol_v_tech_variation[[3]], metric==div_metric & 
                         normalization_type %in% c("none","abundance_based")), 
                  aes(x=variation_label, y=value, fill=variation))+
            geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(div_metric)+
            scale_fill_manual(values = c("#d53e4f", "#3288bd"))
        
        print(p)
        print(q)
        
    } else {
        
        p<-ggplot(subset(biol_v_tech_variation[[3]], metric==div_metric),
                  aes(x=variation_label, y=value, fill=variation))+
            geom_boxplot()+theme_bw()+ theme(axis.text.x=element_text(angle = -45, hjust = 0))+
            facet_grid(pipe~normalization)+xlab("Type of Variation") +ylab(div_metric)+
            scale_fill_manual(values = c("#d53e4f", "#3288bd"))
        
        print(p)
    }
    
    
}

biol_v_tech_variation_summary<-biol_v_tech_variation[[3]] %>% 
    group_by(variation, variation_label, pipe, normalization, metric, normalization_type) %>%
    summarise(mean_value=mean(value), stdev=sd(value), N=length(value), se=stdev/sqrt(N))

print(ggplot(biol_v_tech_variation_summary, 
             aes(x=variation_label, y=normalization, fill=mean_value))+
          geom_tile()+theme_bw()+facet_grid(metric~pipe, scales="free", space="free")+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+
          scale_fill_distiller(palette = "YlGn"))

biol_v_tech_variation_summary$metric<-as.character(biol_v_tech_variation_summary$metric)

biol_v_tech_variation_summary[which(biol_v_tech_variation_summary$metric == "bray" & 
           biol_v_tech_variation_summary$normalization_type=="abundance_based"), c("metric")]<-c("bray_abund")

biol_v_tech_variation_summary[which(biol_v_tech_variation_summary$metric == "wunifrac" & 
                                        biol_v_tech_variation_summary$normalization_type=="abundance_based"), c("metric")]<-c("wunifrac_abund")

biol_v_tech_variation_summary$metric<-factor(biol_v_tech_variation_summary$metric, 
                                             levels=c("jaccard","unifrac", "wunifrac", "bray", "wunifrac_abund", "bray_abund"),
                                             ordered = T)
print(ggplot(biol_v_tech_variation_summary, 
       aes(x=variation_label, y=normalization, fill=mean_value))+
    geom_tile()+theme_bw()+facet_grid(metric~pipe+variation, scales="free", space="free")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))+
    scale_fill_distiller(palette = "YlGn"))

biol_v_tech_variation_summary_highlevel<-biol_v_tech_variation[[3]] %>% 
    group_by(variation, pipe, normalization, metric, normalization_type) %>%
    summarise(mean_value=mean(value), stdev=sd(value), N=length(value), se=stdev/sqrt(N))

print(ggplot(biol_v_tech_variation_summary_highlevel, 
             aes(x=variation, y=normalization, fill=mean_value))+
          geom_tile()+theme_bw()+facet_grid(metric~pipe, scales="free", space="free")+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+
          scale_fill_distiller(palette = "YlGn"))

dev.off()

pdf("./reports/biological_v_technical_variation_diversity_plots_high_level.pdf", height=5, width=8)

print(ggplot(subset(biol_v_tech_variation_summary_highlevel, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=normalization, y=mean_value, color=pipe, linetype=variation, 
                 shape=variation, group=interaction(pipe,variation)))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
            xlab("Rarefaction Level")+ ylab("Mean Distance Value"))

print(ggplot(subset(biol_v_tech_variation_summary_highlevel, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=normalization, y=mean_value, color=pipe, linetype=variation, 
                 shape=variation, group=interaction(pipe,variation)))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se, width=0.5))+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Rarefaction Level")+ ylab("Mean Distance Value"))

print(ggplot(subset(biol_v_tech_variation_summary_highlevel, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=normalization, y=mean_value, color=pipe, linetype=variation, 
                 shape=variation, group=interaction(pipe,variation)))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_wrap(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se, width=0.5))+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Rarefaction Level")+ ylab("Mean Distance Value"))

print(ggplot(subset(biol_v_tech_variation_summary_highlevel, 
                    normalization_type %in% c("none", "abundance_based") &
                        metric %in% c("bray", "wunifrac")), 
             aes(x=interaction(variation, normalization), y=mean_value, color=pipe, 
                 shape=variation, group=interaction(pipe,normalization)))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se, width=0.5))+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Normalization Method")+ ylab("Mean Distance Value"))

print(ggplot(subset(biol_v_tech_variation_summary_highlevel, 
                    normalization_type %in% c("none", "abundance_based") &
                        metric %in% c("bray", "wunifrac")), 
             aes(x=variation, y=mean_value, color=pipe, 
                 shape=variation, group=interaction(pipe,normalization)))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(metric~normalization)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          geom_errorbar(aes(ymin=mean_value-se, ymax=mean_value+se, width=0.5))+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Normalization Method")+ ylab("Mean Distance Value"))

biol_v_tech_variation_summary_higherlevel<- biol_v_tech_variation_summary_highlevel[,1:6] %>% spread(variation, mean_value)

biol_v_tech_variation_summary_higherlevel$mean_difference<-(biol_v_tech_variation_summary_higherlevel$biological - biol_v_tech_variation_summary_higherlevel$technical)

print(ggplot(subset(biol_v_tech_variation_summary_higherlevel, 
                    normalization_type %in% c("none", "abundance_based") &
                        metric %in% c("bray", "wunifrac")), 
             aes(x=normalization, y=mean_difference, color=pipe, group=pipe))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Normalization Method")+ ylab("Mean Biological Distance - Mean Technical Distance"))

 print(ggplot(subset(biol_v_tech_variation_summary_higherlevel, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=normalization, y=mean_difference, color=pipe, group=pipe))+
          geom_point(size=2)+geom_line()+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Rarefaction Level")+ ylab("Mean Biological Distance - Mean Technical Distance"))

biol_v_tech_variation_overall_mean<-biol_v_tech_variation[[3]] %>% 
    group_by(pipe, normalization, metric, normalization_type) %>%
     summarise(overall_mean=mean(value))
 
biol_v_tech_variation_ma_plot<-merge(biol_v_tech_variation_summary_higherlevel, biol_v_tech_variation_overall_mean)

print(ggplot(subset(biol_v_tech_variation_ma_plot, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=overall_mean, y=mean_difference, color=pipe, shape=normalization))+
          geom_point(size=3)+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Overall Mean Distance")+ ylab("Mean Biological Distance - Mean Technical Distance"))

print(ggplot(subset(biol_v_tech_variation_ma_plot, 
                    normalization_type %in% c("none", "rarefaction")), 
             aes(x=overall_mean, y=mean_difference, color=pipe, shape=normalization))+
          geom_point(size=3)+theme_bw()+facet_grid(~metric, scales = "free")+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Overall Mean Distance")+ ylab("Mean Biological Distance - Mean Technical Distance"))

print(ggplot(subset(biol_v_tech_variation_ma_plot, 
                    normalization_type %in% c("none", "abundance_based") &
                        metric %in% c("bray", "wunifrac")), 
             aes(x=overall_mean, y=mean_difference, color=pipe, shape=normalization))+
          geom_point(size=3)+theme_bw()+facet_grid(~metric)+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Overall Mean Distance")+ ylab("Mean Biological Distance - Mean Technical Distance"))

print(ggplot(subset(biol_v_tech_variation_ma_plot, 
                    normalization_type %in% c("none", "abundance_based") &
                        metric %in% c("bray", "wunifrac")), 
             aes(x=overall_mean, y=mean_difference, color=pipe, shape=normalization))+
          geom_point(size=3)+theme_bw()+facet_grid(~metric, scales="free")+
          theme(axis.text.x=element_text(angle = -45, hjust = 0))+ 
          scale_colour_brewer(palette = "Dark2")+ 
          ggtitle("Biological vs. Technical Variation")+
          xlab("Overall Mean Distance")+ ylab("Mean Biological Distance - Mean Technical Distance"))
dev.off()

pdf("./reports/biological_v_technical_variation_diversity_varpart_stats.pdf", height=14, width=16)

varpart_stats$normalization<-factor(varpart_stats$normalization, levels=c("CSS", "TMM", "TSS", "RLE", "UQ", "RAW", "rare10000", 
                                                                          "rare5000", "rare2000", "rareq15"), ordered = T)
varpart_stats$metric<-factor(varpart_stats$metric, levels=c("jaccard", "unifrac", "wunifrac", "bray"), ordered = T)
varpart_stats$txt_label<-paste0(round(varpart_stats$Adj.R.square, 2), "\n", varpart_stats$significance)

print(ggplot(subset(varpart_stats, effect=="conditional"), aes(x=pipe, y=normalization, fill=Adj.R.square))+geom_tile()+facet_grid(metric~feature, scales="free", space="free")+
    theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
    scale_fill_distiller(palette = "YlGn", direction = 1)+ geom_text(aes(label=txt_label))+ 
        ggtitle("Biological vs. Technical Variation varpart Stats: Conditional (partial) Variation"))

print(ggplot(subset(varpart_stats, effect=="marginal"), aes(x=pipe, y=normalization, fill=Adj.R.square))+geom_tile()+facet_grid(metric~feature, scales="free", space="free")+
          theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
          scale_fill_distiller(palette = "YlGn", direction = 1)+ geom_text(aes(label=txt_label))+ 
          ggtitle("Biological vs. Technical Variation varpart Stats: Marginal (simple) Variation"))

print(ggplot(subset(varpart_stats, effect=="global"), aes(x=pipe, y=normalization, fill=Adj.R.square))+geom_tile()+facet_grid(~metric, scales="free", space="free")+
          theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
          scale_fill_distiller(palette = "YlGn", direction = 1)+ geom_text(aes(label=txt_label))+ 
          ggtitle("Biological vs. Technical Variation varpart Stats: Global Variation"))

dev.off()

pdf("./reports/biological_v_technical_variation_diversity_varpart_stats2.pdf", height=8, width=12)

varpart_stats_sub <- subset(varpart_stats, effect %in% c("conditional", "global")) 
varpart_stats_sub$feature<-factor(varpart_stats_sub$feature,
                                      levels=c("all","subject","titration","seq_run"),
                                      ordered=T)
varpart_stats_sub$metric<-factor(varpart_stats_sub$metric,
                                     levels=c("jaccard","unifrac","wunifrac","bray"),
                                     ordered=T)
varpart_stats_sub$Type<-varpart_stats_sub$feature
varpart_stats_sub$Type<-gsub(varpart_stats_sub$Type, pattern="subject", replacement="biological")
varpart_stats_sub$Type<-gsub(varpart_stats_sub$Type, pattern="titration", replacement="biological")
varpart_stats_sub$Type<-gsub(varpart_stats_sub$Type, pattern="seq_run", replacement="technical")
varpart_stats_sub$Type<-gsub(varpart_stats_sub$Type, pattern="all", replacement=" global")

varpart_stats_sub$normalization<-factor(varpart_stats_sub$normalization, 
                                    levels=c("CSS", "TMM", "TSS", 
                                             "RLE", "UQ", "RAW", "rare10000",
                                             "rare5000", "rare2000", "rareq15"), 
                                    ordered = T)

for(value in levels(varpart_stats_sub$metric)){
    print(ggplot(subset(varpart_stats_sub, metric==value), 
                 aes(x=feature, y=Adj.R.square, fill=Type))+
              geom_bar(stat="identity")+facet_grid(pipe~normalization)+
              theme_bw()+theme(axis.text.x=element_text(angle = -45, hjust = 0))+
              ggtitle(paste0("Biological vs. Technical Variation varpart Stats: ", value))+
              xlab("")+ylab("Adjusted R^2")+
              scale_fill_manual(values = c("gray","#d53e4f", "#3288bd")))
}

dev.off()
