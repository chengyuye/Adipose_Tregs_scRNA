library(readxl)
library(tidyverse)
library(patchwork)
library(viridis)
library(Seurat)
library(scCustomize)
library(ggsci)
library(RColorBrewer)
## density plots function
library(Nebulosa)
library(stringr)

## set up working directory
setwd("/Users/chengyuye/Library/CloudStorage/OneDrive-EmoryUniversity/Rotation_2_Li Lab/scRNA_TCR/Data/processed_data/")

## read in data
adipose.harmony.subset <- readRDS('../processed_data/adipose/seurat_update/subset_adipose_harmony_final_res0.4_mostrecent.rds')

ch.markers <- read_excel('../processed_data/adipose/Main Cholesterol Biosynthesis Branch_Gene List.xlsx')
ch.markers <- ch.markers$`Main Cholesterol Biosynthesis Branch Genes`

ch_main.markers <- read_excel('../processed_data/Core Cholesterol gene signature_Mouse.xlsx',col_names = FALSE)
ch_main.markers <- ch_main.markers$...1

## TCR upregulation signatures
tcr_up <- read_excel('../processed_data/adipose/TCR_activation_upregulated_genes.xlsx')
tcr_up <- tcr_up$Uhrf1


## TCR downregulation signatures
tcr_down <- read_excel('../processed_data/adipose/TCR_activation_downregulated_genes.xlsx',col_names = FALSE)
tcr_down <- tcr_down$...1

## VAT Treg up-regulated genes
vat_up <- read.table('./adipose/Fat Treg specific signature up gene list-DC.txt')
vat_up <- vat_up$V1

## VAT Treg down-regulated genes
vat_down <- read.table('./adipose/Fat Treg specific signature down gene list-DC.txt')
vat_down <- vat_down$V1


## IFNa signatures
ifna_up <- read.table('./adipose/IFNa signature up-signature.txt')
ifna_up <- ifna_up$V1
ifna_up <- str_to_title(ifna_up)


## IFNg signatures
ifng_up <- read.table('./adipose/IFNg signature up-signature.txt')
ifng_up <- ifng_up$V1
ifng_up <- str_to_title(ifng_up)
ifng_up <- append(ifng_up, 'H2-DMa')


## CD44 hi WT vs KO signatures
cd44_wt_ko <- read.table('./adipose/CD44hi WT vs CD44hi KO_up signature_gene list.txt')
cd44_wt_ko <- cd44_wt_ko$V1



adipose.harmony.subset$celltype_group <- paste(adipose.harmony.subset$seurat_clusters,
                                          sep = '_',
                                          adipose.harmony.subset$Group)

### addmodulescore
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(ch.markers),
                                         name = 'Cholesterol_signature')

adipose.harmony.subset$Group <- as.factor(adipose.harmony.subset$Group)


## main markers module
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(ch_main.markers),
                                         name = 'Main_Cholesterol_signature')

## main markers module
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(tcr_up),
                                         name = 'TCR_activation_signature')


## main markers module
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(tcr_down),
                                         name = 'TCR_down_signature')


## VAT_Treg up 
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(vat_up),
                                         name = 'VAT_upregulated')


## VAT_Treg down
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(vat_down),
                                         name = 'VAT_downregulated')


## VAT_ifna up
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(ifna_up),
                                         name = 'Ifna_upregulated')

## VAT_ifng up
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(ifng_up),
                                         name = 'Ifng_upregulated')

## CD44 up
adipose.harmony.subset <- AddModuleScore(adipose.harmony.subset,
                                         features = list(cd44_wt_ko),
                                         name = 'CD44_upregulated')


## another color scheme to consider
col <- scale_color_gradientn(colors = c("grey85", brewer.pal(7, "Reds")))
###for dimplot
colorRampPalette(brewer.pal(n = 10, name = "Paired"))(5)

#1. not split by condition
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'Cholesterol_signature1', label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)
# save the plot
ggsave('Cholesterol_pathway_overall_umap_updated.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 6, height = 5) 



#2.  split by condition
p <- FeaturePlot(adipose.harmony.subset,
                 split.by = 'Group',
                 features = 'Cholesterol_signature1', label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)
# save the plot
ggsave('Cholesterol_pathway_split_by_groups_umap_updated.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 5) 


#3.  split by condition ------ Main cholesterol signatures
p <- FeaturePlot(adipose.harmony.subset,
                 split.by = 'Group',
                 features = 'Main_Cholesterol_signature1', label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)
# save the plot
ggsave('Main_Cholesterol_pathway_split_by_groups_umap_updated.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 5) 


#4.  not spliting ------ Main cholesterol signatures
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'Main_Cholesterol_signature1', label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)
# save the plot
ggsave('Main_Cholesterol_pathway_umap_updated1.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 6, height = 5) 


#5.  not spliting ------ TCR activation signatures
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'TCR_activation_signature1', label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)

# save the plot
ggsave('TCR_activation_umap.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 6, height = 5) 

#6.  spliting ------ TCR activation signatures
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'TCR_activation_signature1', 
                 split.by = 'Group',
                 label = TRUE, repel = TRUE,
                 pt.size = 1,label.size = 5,
                 cols=c("grey","yellow","red","brown")
)

# save the plot
ggsave('TCR_activation_umap_split_by_group.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 10, height = 5) 

### using scCustom
FeaturePlot_scCustom(seurat_object = adipose.harmony.subset, 
                    features = "Cholesterol_signature1", split.by = 'Group',na_cutoff = NA)


### scCustom Stacked_VlnPlot
p <- VlnPlot(adipose.harmony.subset, 
        features = "Cholesterol_signature1",
        split.by = 'Group')
ggsave('Cholesterol_pathway_split_by_groups_violinplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 5) 

p <- VlnPlot(adipose.harmony.subset, 
             features = "Srebf2",
             split.by = 'Group')
ggsave('Srebf2_split_by_groups_violinplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 5) 

sample_colors <- c("dodgerblue",  "firebrick1")
Stacked_VlnPlot(seurat_object = adipose.harmony.subset, 
                features = "Cholesterol_signature1", 
                x_lab_rotate = TRUE,
                colors_use = sample_colors, 
                split.by = "Group")
## colors 
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
                        '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
                        '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398',
                        '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B',
                        '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
                        '#968175'
)

# rename the cell to avoid errors
#colnames(mtx$Ara) <- paste("Ara", colnames(mtx$Ara), sep = "_")
adipose.harmony.subset$celltype_group <- paste(adipose.harmony.subset$seurat_clusters, 
                                               adipose.harmony.subset$Group, sep = "_")


#### dotplot ---- cholesterol signatures
p <- DotPlot(adipose.harmony.subset,
        group.by = 'celltype_group',
        features = 'Cholesterol_signature1',
        dot.scale = 8)+
  scale_color_viridis()
        #cols = my36colors)
ggsave('Cholesterol_pathway_split_by_groups_dotplot1.pdf', p, path = "../../Figures/adipose_updated/",
       width = 7, height = 6) 

## dotplot ---- main cholesterol signatures
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Main_Cholesterol_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('Main_Cholesterol_pathway_split_by_groups_dotplot1.pdf', p, path = "../../Figures/adipose_updated/",
       width = 7, height = 6) 

p <- DotPlot(adipose.harmony.subset,
             features = 'Main_Cholesterol_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('Main_Cholesterol_pathway_dotplot1.pdf', p, path = "../../Figures/adipose_updated/",
       width = 7, height = 6)

### TCR activation signatures ----- upregulation
## dotplot 
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'TCR_activation_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('TCR_activation_dotplot_split_by_group.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6) 

p <- DotPlot(adipose.harmony.subset,
             features = 'TCR_activation_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('TCR_activation_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6)


### TCR activation signatures ----- downregulation
## dotplot 
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'TCR_down_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('TCR_down_dotplot_split_by_group.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6) 

p <- DotPlot(adipose.harmony.subset,
             features = 'TCR_down_signature1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('TCR_down_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6)



### CD44 WT vs KO up signatures
## dotplot 
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'CD44_upregulated1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('CD44_wt_ko_dotplot_split_by_group.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6) 

p <- DotPlot(adipose.harmony.subset,
             features = 'CD44_upregulated1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('CD44_wt_ko_activation_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6)



## dotplot -----srebf2
p <- DotPlot(adipose.harmony.subset,
             features = 'Srebf2',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)

## dotplot -----Pparg
p <- DotPlot(adipose.harmony.subset,
             features = 'Pparg',
             group.by = 'celltype_group',
             dot.scale = 8)+
  scale_color_viridis()

ggsave('Pparg_split_by_group_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 

## dotplot -----Nr4a1
p <- DotPlot(adipose.harmony.subset,
             features = 'Nr4a1',
             group.by = 'celltype_group',
             dot.scale = 8)+
  scale_color_viridis()

ggsave('Nr4a1_split_by_group_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 6, height = 6) 

### Il27ra ---- overall 
p <- DotPlot(adipose.harmony.subset,
             features = 'Il27ra',
             dot.scale = 8)+
  scale_color_viridis()

ggsave('Il27ra_overall_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 

### Il27ra ---- split by group
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Il27ra',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('Il27ra_split_by_groups_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 

### Il27ra  ---- overall ----umap
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'Il27ra',
                 pt.size = 0.6)

ggsave('Il27ra_overall_featureplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 7, height = 6) 

### Il27ra ---- split by group---umap
p <- FeaturePlot(adipose.harmony.subset,
                 features = 'Il27ra',
                 pt.size = 0.6,
                 split.by = 'Group')
#cols = my36colors)
ggsave('Il27ra_split_by_groups_featureplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 10, height = 5) 

### Ifngr1  ---- overall 
p <- DotPlot(adipose.harmony.subset,
             features = 'Ifngr1',
             dot.scale = 8)+
  scale_color_viridis()

ggsave('Ifngr1_overall_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 

### Ifngr1 ---- split by group
p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Ifngr1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('Ifngr1_split_by_groups_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 

### Ifngr1  ---- overall ----umap
p <- FeaturePlot(adipose.harmony.subset,
              features = 'Ifngr1',
             pt.size = 0.6)

ggsave('Ifngr1_overall_featureplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 7, height = 6) 

### Ifngr1 ---- split by group---umap
p <- FeaturePlot(adipose.harmony.subset,
             features = 'Ifngr1',
             pt.size = 0.6,
             split.by = 'Group')
#cols = my36colors)
ggsave('Ifngr1_split_by_groups_featureplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 10, height = 5) 





p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Srebf1',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)


p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Fdft1',
             dot.scale = 8)+
  scale_color_viridis()


p <- DotPlot(adipose.harmony.subset,
             group.by = 'celltype_group',
             features = 'Srebf2',
             dot.scale = 8)+
  scale_color_viridis()
#cols = my36colors)
ggsave('Srebf2_split_by_groups_dotplot.pdf', p, path = "../../Figures/adipose_updated/",
       width = 6, height = 6) 


### Dotplot for all the genes
dot_plot <- function(marker){
  p <- DotPlot(adipose.harmony.subset, 
               group.by  = 'celltype_group',
               features = as.character(marker),
               dot.scale = 8)+
       scale_color_viridis()
  ggsave(paste0(as.character(marker), '_split_by_group_dotplot.pdf'), p,
         path = "../../Figures/adipose_updated/cholesterol pathway/",
         width = 6, height = 6)
}

lapply(ch.markers, dot_plot)




### density plots for the VAT Treg upregulated genes
## overall
p <- plot_density(adipose.harmony.subset, 
             features = 'VAT_upregulated1', 
             reduction = 'umap',
             method = 'wkde')

ggsave('VAT_Treg_upregulated_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7.25, height = 6) 

## split by group
p <- plot_density(adipose.harmony.subset, 
                  features = 'VAT_upregulated1', 
                  reduction = 'umap',
                  method = 'wkde') + 
  facet_grid(.~adipose.harmony.subset$Group)

ggsave('VAT_Treg_upregulated_splitbygroup_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 6) 


### density plots for the VAT Treg downregulated genes
## overall
p <- plot_density(adipose.harmony.subset, 
                  features = 'VAT_downregulated1', 
                  reduction = 'umap',
                  method = 'wkde')

ggsave('VAT_Treg_downregulated_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7.25, height = 6) 

## split by group
p <- plot_density(adipose.harmony.subset, 
                  features = 'VAT_downregulated1', 
                  reduction = 'umap',
                  method = 'wkde')+ 
  facet_grid(.~adipose.harmony.subset$Group)

ggsave('VAT_Treg_downregulated_splitbygroup_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 6) 

### density plots for the Ifna Treg upregulated genes
## overall
p <- plot_density(adipose.harmony.subset, 
                  features = 'Ifna_upregulated1', 
                  reduction = 'umap',
                  method = 'wkde')

ggsave('VAT_Ifna_upregulated_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7.25, height = 6) 

## split by group
p <- plot_density(adipose.harmony.subset, 
                  features = 'Ifna_upregulated1', 
                  reduction = 'umap',
                  method = 'wkde')+ 
  facet_grid(.~adipose.harmony.subset$Group)

ggsave('VAT_Ifna_upregulated_splitbygroup_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 6) 


### density plots for the Ifng Treg upregulated genes
## overall
p <- plot_density(adipose.harmony.subset, 
                  features = 'Ifng_upregulated1', 
                  reduction = 'umap',
                  method = 'wkde')

ggsave('VAT_Ifng_upregulated_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7.25, height = 6) 

## split by group
p <- plot_density(adipose.harmony.subset, 
                  features = 'Ifng_upregulated1', 
                  reduction = 'umap',
                  method = 'wkde')+ 
  facet_grid(.~adipose.harmony.subset$Group)

ggsave('VAT_Ifng_upregulated_splitbygroup_densityplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 12, height = 6) 



### dot plots for the VAT Treg upregulated genes
## split-by group
p <- DotPlot(adipose.harmony.subset, 
                  features = 'VAT_upregulated1', 
                  group.by = 'celltype_group')+scale_color_viridis()
 
ggsave('VAT_Treg_upregulated_splitbygroup_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7.25, height = 6) 






### dot plots for the Ifna upregulated genes
## split by group
p <- DotPlot(adipose.harmony.subset, 
                  features = 'Ifna_upregulated1', 
             group.by = 'celltype_group')+scale_color_viridis()

ggsave('Ifna_upregulated_splitbygroup_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6) 




### dot plots for the Ifng upregulated genes
## split by group
p <- DotPlot(adipose.harmony.subset, 
             features = 'Ifng_upregulated1', 
             group.by = 'celltype_group')+scale_color_viridis()

ggsave('Ifng_upregulated_splitbygroup_dotplot.pdf', p, path = "../../Figures/adipose_updated/cholesterol pathway/",
       width = 7, height = 6) 

