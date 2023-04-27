library(tidyverse)
library(SeuratDisk)
library(Seurat)
library(anndata)


options(scipen=999)
options(repr.matrix.max.rows=20)
options(repr.matrix.max.cols=50)


library(cowplot)

library(patchwork)
library(ggsci)


library(viridis)
library(scales)


library(ggpubr)
library(scales)


text_sizes <- theme(axis.text.x=element_text(size=9,colour="black"),
                      axis.text.y=element_text(size=9,colour="black"),
                      axis.title.y=element_text(size=9,colour="black", margin = margin(t = 2, l = 1, r = 1, b =2, unit = "pt")),
                      axis.title.x=element_text(size=9,colour="black", margin = margin(t = 1, l = 2, r = 2, b = 1, unit = "pt")),
                      legend.text = element_text(size=9,colour="black"),
                      legend.title = element_text(size=9,colour="black", margin = margin(t = 5, l = 0, r = 0, b = 5, unit = "pt")),
                      legend.key = element_rect(colour="transparent", fill = "transparent"),
                      strip.text.x = element_text(size=6,color = 'black',face="bold", angle=0),
                      strip.text.y = element_text(size=6,color = 'black', face="bold", angle=0, vjust=0.5, hjust=0),
                      axis.ticks= element_line(color = 'black', size=0.2),
                      axis.line = element_line(colour = "gray50", size = 0.2, linetype = "solid"),
                      plot.margin=unit(c(0,0,0,0),"pt"),
                      plot.title=element_text(size=7, face='plain', colour="black"))


common_minimal <- text_sizes + theme(
        plot.background = element_rect(fill = NA,colour = NA),
        strip.background = element_rect(fill = NA,colour = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) + theme(
    legend.spacing = unit(0.15, 'cm'), 
    legend.key.size = unit(0.2, "cm"))


# commonly used, x axis text 45 degree
common_0x <- common_minimal + theme(axis.text.x = element_text(angle=0))

color_pal = readRDS('/nfs/research/irene/ysong/heart_cell_type_color_pal_5.rds')
color_pal_species = readRDS('/nfs/research/irene/ysong/heart_species_color_pal_5.rds')


metrics = readr::read_csv("/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/overall_rank_figures_noLISI_noSAMap/Heart_hs_mf_overall_rank.csv")

metrics <- metrics %>% filter(integration_method_batch != 'unintegrated')
metrics <- metrics %>% arrange(type)
myplots <- list()


for(i in seq(1, nrow(metrics))){
    
type_now = metrics$type[i]
method_new = metrics$integration_method_batch[i]

    message(type_now)
if(grepl("SH", type_now)) {
    
    homo_name = 'many_higher_homology_conf'
} else if (grepl("HE", type_now)) {
    
    homo_name = 'many_higher_expr'
} else if (grepl("O2O", type_now)) {
    
    homo_name = 'one2one_only'
}

    
method_old = method_new
    
method_old = gsub("SeuratV4 CCA", "seuratCCA", method_old)
method_old = gsub( "SeuratV4 RPCA", "seuratRPCA", method_old)
method_old = gsub( "Scanorama", "scanorama", method_old)
method_old = gsub( "LIGER UINMF","rligerUINMF", method_old)
method_old = gsub("Harmony", "harmony", method_old)


SM_score = metrics %>% filter(type == type_now) %>% pull(avg_score_batch) %>% round(2)
BC_score = metrics %>% filter(type == type_now) %>% pull(avg_score_bio) %>% round(2)
    
if(method_old != 'rligerUINMF'){
    
    file = paste0("/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/Heart_hs_mf/results/", method_old, "/cross_species/integrated_h5ad/metadata_nf_", homo_name, "_", method_old, "_integrated.h5seurat")
    
    
} else {
    
    file = paste0("/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/Heart_hs_mf/results/", method_old, "/cross_species/integrated_h5ad/rliger_uinmf_metadata_", homo_name, "_", method_old, "_integrated.h5seurat")
    
}
    
message(file.exists(file))
    


obj <- LoadH5Seurat(file, assays = 'RNA',reductions = 'umap',misc = FALSE, verbose = FALSE)
                    

cell_type_pal =color_pal

species_pal = hue_pal()(3)
names(species_pal) = levels(factor(obj@meta.data$species))

species = DimPlot(object = obj, reduction = 'umap', group.by = c('species'), pt.size = 0.02, shuffle = TRUE) + 
theme(legend.position = 'none') + text_sizes + labs(title = paste0(type_now, "\n\n", "SM: ", SM_score, " BC: ", BC_score)) + scale_color_manual(values = species_pal) +
theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

cell_type = DimPlot(object = obj, reduction = 'umap', group.by = 'cell_ontology_mapped', pt.size = 0.02, shuffle = TRUE) +
theme(legend.position = 'none') + text_sizes + theme(plot.title = element_blank()) + scale_color_manual(values = cell_type_pal)+
theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

fig = plot_grid(species, cell_type, ncol = 1, rel_heights  = c(1, 0.9))
                    
#png(paste0("panc_", type_now, "_umap.png"), height = 4, width = 2, res = 300, unit = 'in')
#print(fig)
#dev.off()
    
myplots[[i]] <- fig

}


samap <- "/nfs/research/irene/ysong/RESULTS/NEXTFLOW/CROSS-SPECIES_test/Heart_hs_mf/results/SAMap/cross_species/integrated_h5ad/metadata_nf_many_higher_expr_SAMap_integrated.h5seurat"

type_now = 'SAMap all genes'

obj <- LoadH5Seurat(samap, assays = 'RNA',reductions = 'umap',misc = FALSE, verbose = FALSE)



species = DimPlot(object = obj, reduction = 'umap', group.by = c('species'), pt.size = 0.02, shuffle = TRUE) + 
theme(legend.position = 'none') + text_sizes + labs(title = paste0(type_now, "\n\n\n")) + scale_color_manual(values = species_pal) +
theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

cell_type = DimPlot(object = obj, reduction = 'umap', group.by = 'cell_ontology_mapped', pt.size = 0.02, shuffle = TRUE) +
theme(legend.position = 'none') + text_sizes + theme(plot.title = element_blank()) + scale_color_manual(values = cell_type_pal)+
theme(axis.line = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank())

fig = plot_grid(species, cell_type, ncol = 1, rel_heights  = c(1, 0.9))
                    
#png(paste0("panc_", type_now, "_umap.png"), height = 4, width = 2, res = 300, unit = 'in')
#print(fig)
#dev.off()
message('samap')
    
myplots[[i+1]] <- fig


saveRDS(myplots, 'heart_hs_mf_all_plots_samap.rds')


png(paste0("heart_hs_mf_umap_all_one_samap.png"), height = 8, width = 7.5, res = 300, unit = 'in')

print(plot_grid(myplots[1][[1]], myplots[2][[1]], myplots[3][[1]], myplots[4][[1]], myplots[5][[1]], myplots[6][[1]], 
                myplots[7][[1]], myplots[8][[1]], myplots[9][[1]], myplots[10][[1]], myplots[11][[1]], myplots[12][[1]]
                , myplots[13][[1]], myplots[14][[1]], myplots[15][[1]], myplots[16][[1]]
                , myplots[17][[1]], myplots[18][[1]], ncol = 6, axis = 'lb'))

dev.off()

png(paste0("heart_hs_mf_umap_all_two_samap.png"), height = 6, width = 7.5, res = 300, unit = 'in')

print(plot_grid(myplots[19][[1]], myplots[20][[1]], myplots[21][[1]], myplots[22][[1]],myplots[23][[1]], myplots[24][[1]], myplots[25][[1]],myplots[26][[1]], myplots[27][[1]], myplots[28][[1]],ncol = 6, axis = 'lb'))

dev.off()


png(paste0("heart_hs_mf_umap_all_both_samap.png"), height = 14, width = 7.5, res = 300, unit = 'in')

print(plot_grid(myplots[1][[1]], myplots[2][[1]], myplots[3][[1]], myplots[4][[1]], myplots[5][[1]], myplots[6][[1]], 
                myplots[7][[1]], myplots[8][[1]], myplots[9][[1]], myplots[10][[1]], myplots[11][[1]], myplots[12][[1]],
                myplots[13][[1]], myplots[14][[1]], myplots[15][[1]], myplots[16][[1]], myplots[17][[1]], myplots[18][[1]],
                myplots[19][[1]], myplots[20][[1]], myplots[21][[1]], myplots[22][[1]],myplots[23][[1]], myplots[24][[1]], 
                myplots[25][[1]],myplots[26][[1]], myplots[27][[1]], myplots[28][[1]],ncol = 6, axis = 'lb'))

dev.off()