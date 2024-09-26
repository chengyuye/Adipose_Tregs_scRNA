library(Seurat)
library(harmony)
library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(umap)
library(reticulate)
library(scRepertoire)
library(ggsci)
library(ggpubr)

theme_bar <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),#去网格线
          axis.title.x = element_blank(),#去x轴标签
          axis.title.y=element_text(face = "bold",size = 14),#y轴标签加粗及字体大小
          axis.text = element_text(face = "bold",size = 12),#坐标轴刻度标签加粗
          # axis.ticks = element_line(color='black'),#坐标轴刻度线
          # axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),#去除图例标题
          # legend.justification=c(1,0),#图例在画布的位置(绘图区域外)
          legend.position=c(0.75, 0.9),#图例在绘图区域的位置
          # legend.position='top',#图例放在顶部
          legend.direction = "horizontal",#设置图例水平放置
          # legend.spacing.x = unit(2, 'cm'),
          legend.text = element_text(face = "bold",size = 12,margin = margin(r=20)),
          legend.background = element_rect( linetype="solid",colour ="black")
          # legend.margin=margin(0,0,-7,0)#图例与绘图区域边缘的距离
          # legend.box.margin =margin(-10,0,0,0)
    )
  
}
## check the cell number
table(adipose.harmony.subset$Group)

## check the cell number per cluster 
table(Idents(adipose.harmony.subset), adipose.harmony.subset$Group)#


prop.table(table(Idents(adipose.harmony.subset)))

Cellratio <- prop.table(table(Idents(adipose.harmony.subset), adipose.harmony.subset$Group), margin = 2)#
Cellratio


Cellratio <- as.data.frame(Cellratio)
colnames(Cellratio)[1] <- 'celltype'
colnames(Cellratio)[2] <- 'sample'

## save the cell ratio result
write.csv(Cellratio, file = "C:/BGI/mianhua/CK_harmony_updated_test/dim50_annotation/DA/cell_ratio.csv")


library(ggplot2)
p <- ggplot(data = Cellratio,aes(x =celltype, y= Freq, fill = sample)) + 
  geom_bar(stat = "identity",position = "dodge2", width = 0.8)+ 
  theme_classic() +
  labs(x='Celltype',y = 'Ratio') +
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_bar()
#coord_flip()+
#theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p <- p +  scale_fill_npg()
ggsave('cellratio_barplot_split_celltype.pdf', p, path = "../../Figures/adipose_updated/DA/",
       width = 8, height = 6)


ggbarplot(Cellratio, x = "celltypes", y = "Freq", add = "mean_se",
          color = "supp", palette = "jco", 
          position = position_dodge(0.8))+
  stat_compare_means(aes(group = supp), label = "p.signif", label.y = 29)


p <- ggplot(Cellratio) + 
    geom_bar(aes(x =sample, y= Freq, fill = celltype),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+ 
    theme_classic() +
    labs(x='Sample',y = 'Ratio')+
    #scale_fill_manual(values = allcolour)+
    theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))+
    scale_fill_d3()


ggsave('cellratio_barplot_celltype.pdf', p, path = "../../Figures/adipose_updated/DA/",
       width = 6, height = 6)
