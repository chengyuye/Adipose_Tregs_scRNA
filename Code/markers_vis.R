library(Seurat)
library(harmony)
library(future)
library(reticulate)
library(umap)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(patchwork)
library(data.table)
library(readxl)
library(monocle)
library(ggsci)
library(ggpubr)
library(patchwork)
library(clusterProfiler)
library(MySeuratWrappers)  



markers <- c('Tbx21', 'Cxcr3', 'Tcf7', 'Nt5e', 
             'Dusp10', 'Icos', 'Pparg', 'Il1rl1', 
             'Areg', 'Nr4a1', 'Ly6c1', 'Txnip',
             'Cd74')

markers <- c('Il1rl1', 'Nt5e', 'Tcf7', 'Tbx21')


#??ɫ????UMAPͼһ??
col <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175')#??ɫ????  
p <- VlnPlot(adipose.harmony.subset, features = markers,  
             stack=T, pt.size=0,  
             #cols = col,#??ɫ  
             direction = "horizontal", #ˮƽ??ͼ
             x.lab = '', y.lab = '')+ 
  theme(axis.text.x = element_blank(),   
        axis.ticks.x = element_blank())#????ʾ?????̶?
p <- p + scale_fill_d3()

ggsave('marker_stacked_violinplot.pdf', p, path = "../../Figures/adipose_updated/marker_final/",
       width = 8, height = 5)

ggsave('4marker_stacked_violinplot.pdf', p, path = "../../Figures/adipose_updated/marker_final/",
       width = 6, height = 4)



###dot plot
p <- DotPlot(adipose.harmony.subset, features = markers, 
             dot.scale = 15) +
  coord_flip() +
  theme_bw()+ #  
  theme(panel.grid = element_blank(),axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#330066','#336699','#66CC66','#FFCC33'))+ 
  labs(x=NULL,y=NULL)+
  guides(size=guide_legend(order=3)) 

ggsave('markers_dotplots.pdf', p, path = "../../Figures/adipose_updated/marker_final/",
       width = 6, height = 8)


### feature plot
feature_plot <- function(marker){
  p <- FeaturePlot(adipose.harmony.subset, features = as.character(marker),pt.size = 0.8)
  ggsave(paste0(as.character(marker), '_featureplot.pdf'), p, path = "../../Figures/adipose_updated/marker_final/",
         width = 7, height = 6)
}

FeaturePlot(adipose.harmony.subset, "Cd74")  
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

lapply(markers, feature_plot)
