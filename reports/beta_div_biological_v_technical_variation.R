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
    summarise(mean_value=mean(value))

print(ggplot(biol_v_tech_variation_summary, 
       aes(x=variation_label, y=normalization, fill=mean_value))+
    geom_tile()+theme_bw()+facet_grid(metric~pipe, scales="free", space="free")+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))+
    scale_fill_distiller(palette = "RdBu"))

dev.off()