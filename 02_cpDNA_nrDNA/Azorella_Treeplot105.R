setwd("C:/Users/7827x/OneDrive/Desktop/Azorella_NewBatch/01_Chapter2/02_GenomeSkimming")

library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ggplot2)
library(aplot)
library(ggpubr)
library(ggbreak)

tax <- read.table("ColoringGroup.csv", sep = ',', header  = T)
colnames(tax) = c("treename", "group", "ploidy")
tax1 <- as.data.frame(tax)

plastome <- read.tree(file = "Azorella_plastome_n99complete_n6map_aln_trim2_rerooted.treefile")

tree1 = full_join(as_tibble(plastome), tax1, by = c('label' = 'treename'))
tree1$treename <- paste0( " [", tree1$ploidy, "]" , tree1$label)
tree2 = as.treedata(tree1)



###############################

ggtree(tree2) +   
  geom_tiplab(aes( label = treename), align = T, show.legend=FALSE) + 
  scale_x_break(c(0.0005, 0.006)) +  hexpand(0.05)

tree1.test <- data.frame(tree1)
tree1.test$bsvalue <- as.numeric(tree1.test$label)

tree1.nodes <- subset(tree1.test, bsvalue >= 90)$node 
       
ggtree(tree2) + 
  geom_text(aes(label=node), hjust=-.3) +
  geom_point2(aes(subset=(node %in% tree1.nodes)),color="green",size=1)


##########################################################

T0 = flip(ggtree(tree2), 29, 28) +      
  geom_tiplab(aes(color = group, label = treename), align = T)  + #, fontface = 4
  geom_point2(aes(subset=(node %in% tree1.nodes)),color="red",size=1.5, alpha = 1) +
  geom_rootedge(rootedge = 0.0005) +
  
  theme(legend.position = "none") + 
  scale_color_manual(values=c("NZ1"= "#D55E00",
                              "NZ2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "SA" = "black")) +
  
  geom_strip("Azcolensoi_MauB1", "Azcockaynei_Man5", "NZ1", color = "#D55E00", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 5, parse = TRUE) +
  
  geom_strip("Azcyanopetala_TeA10", "Azcyanopetala_Tak1", "NZ2_cp1", color = "#0072B2", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azfragosea_CANB797854", "Azfragosea_CANB797887", "Au", color = "#CC79A7", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azcyanopetala_Bro3", "Azsp_AN58", "NZ2_cp2", color = "#0072B2", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025,  
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azschizeilema_Ome", "Azschizeilema_Ome", "Sub", color = "#009E73", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  

  geom_strip("Azpolaris_End2", "Azlyallii_CHR542359", "SubM", color = "#009E73", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azburkartii_NYBG2434", "Azranunculus_NYBG2447", "SA", color = "#000000", 
             offset = 0.0075, align = TRUE, offset.text = 0.00025, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  ggtitle("a) Genome skimming of plastome") +
  geom_treescale(x=0, y =-0.5, fontsize=6, linesize=2, offset=-1.5) +
  scale_y_continuous(expand=c(0, 1.5)) +
  theme(plot.title = element_text(hjust = 0.5, size = 18, vjust = -1)) +
  scale_x_break(c(0.0005, 0.006)) +  hexpand(0.1)


T1 <- T0 +   annotate("segment", x = 0.02057, xend = 0.02057, y = 7.5, yend = 8.5, size =2, colour = "#009E73")

T1

#group <- melt(as.matrix(read.csv("ColoringGroup.csv", header = T, row.names=1)))

######################################################################################
nrDNA <- read.tree(file = "Azorella_nrDNA_sorted_aln_trimgt07_rerooted.treefile")

tree3 = full_join(as_tibble(nrDNA), tax1, by = c('label' = 'treename'))
tree3$treename <- paste0(  tree3$label, " [", tree3$ploidy, "]" )
tree4 = as.treedata(tree3)


###############################

ggtree(tree4) + scale_x_break(c(0.001, 0.02)) +  hexpand(0.05)

tree3.test <- data.frame(tree3)
tree3.test$bsvalue <- as.numeric(tree3.test$label)

tree3.nodes <- subset(tree3.test, bsvalue >= 90)$node 

ggtree(tree4) + 
  scale_x_break(c(0.001, 0.02)) +  hexpand(0.05) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_point2(aes(subset=(node %in% tree3.nodes)),color="green",size=1)


##########################################################

T2 <- ggtree(tree4) +   
  geom_tiplab(hjust =1, aes(color = group, label = treename), align = T) + #fontface = 4
  geom_point2(aes(subset=(node %in% tree3.nodes)),color="red",size=1.5, alpha = 0.9) +
  geom_rootedge(rootedge = 0.005) +
  
  theme(legend.position = "none")  +
#  ggtitle("b) nrDNA Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, vjust = -1)) + 
  scale_color_manual(values=c("NZ1"= "#D55E00",
                              "NZ2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "SA" = "black")) +

  geom_strip("Azcyanopetala_Lux5", "Azhaastii_AN63", "NZ2", color = "#0072B2", 
             offset = 0.047, align = TRUE, offset.text = 0.01, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azallanii_Rau1", "Azroughii_AN67", "NZ1_nr1", color = "#D55E00", 
             offset = 0.047, align = TRUE, offset.text = 0.017, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azschizeilema_Ome", "Azschizeilema_Ome", "Sub", color = "#009E73", 
             offset = 0.047, align = TRUE, offset.text = 0.01, 
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  
  geom_strip("Azfragosea_CANB797887", "Azfragosea_CANB797854", "Au", color = "#CC79A7", 
             offset = 0.047, align = TRUE, offset.text = 0.009,
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azpolaris_Camp2", "Azlyallii_CHR542359", "SubM", color = "#009E73", 
             offset = 0.047, align = TRUE, offset.text = 0.012,
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azcolensoi_MauA5", "Aznitens_AN49", "NZ1_nr2", color = "#D55E00", 
             offset = 0.047, align = TRUE, offset.text = 0.017,
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  
  geom_strip("Azburkartii_NYBG714", "Azranunculus_NYBG2447", "SA", color = "#000000", 
             offset = 0.047, align = TRUE, offset.text = 0.009,
             barsize = 2, fontsize = 4.5, parse = TRUE) +
  ggtitle("b) Genome skimming of nrDNA") +
  scale_y_continuous(expand=c(0, 1.5)) +
  geom_treescale(x=0, y =-0.5, fontsize=6, linesize=2, offset=-1.5) +
  scale_x_break(c(0.001, 0.02)) +  hexpand(0.05) +
  scale_x_reverse()  

T3 <- T2 +   annotate("segment", x = 0.1186, xend = 0.1186, y = 30.5, yend = 31.3, size =2, colour = "#009E73")
T3

png("Fig4_Azorella_plastome_nrDNA3.png", width = 4100, height = 4500, res =250)
ggarrange(T1,T3 #, #label.x = 0.2,
#          labels = c("a) Genome skimming of Plastome", 
#                     "b) Genome skimming of nrDNA"), 
#          hjust = -0.5
)
dev.off()

fig4 <- plot_list(T1, T2)

pdf("Fig4_Azorella_plastome_nrDNA3.pdf", width = 16.8, height = 18)
#fig4
ggarrange(T1,T3#, label.x = 0.2,
          #labels = c("a) Genome skimming of Plastome", 
          #           "b) Genome skimming of nrDNA"), 
         # hjust = -0.5
         )
dev.off()

#scale_color_manual(values=c("Azal_rou"= "#D55E00",
#                            "Azcoc_nit" = "blue", 
#                            "Azcol_ho" = "green",
#                            "Azcya" = "skyblue", 
#                            "Azex" = "antiquewhite", 
#                            "Azfra" = "brown2",
#                            "Azha" = "cornflowerblue",
#                            "Azhy" = "#999999",
#                            "Azlya_po" = "darksalmon",
#                            "Azpa" = "darkseagreen2",
#                            "Azran_bu_lyc" = "deeppink2",
#                            "Azsp" = "black")) 

######################################################################################
######################################################################################
######################################################################################
cladeplastome <- tree_subset(tree2, node=146, levels_back=0)
cladeplastome.bs <- data.frame(as_tibble(cladeplastome))
cladeplastome.bs$bsvalue <- as.numeric(cladeplastome.bs$label)
cladeplastome.bs.node <- subset(cladeplastome.bs, bsvalue >= 90)$node 


T3 <- ggtree(cladeplastome) + 
  geom_point2(aes(subset=(node %in% cladeplastome.bs.node)),color="red",size=2) +
  geom_tiplab(aes(color = group.y), align = T, show.legend = F) +
  scale_color_manual(values=c("G1"= "#D55E00",
                              "G2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "Outgroup" = "black")) +
  scale_y_continuous(expand=c(0, 1.5)) +
  geom_treescale(x=0, y =-0.5, fontsize=6, linesize=2, offset=-1) +
  hexpand(1) 

######################################################################################
######################################################################################
######################################################################################
cladenrDNA <- tree_subset(tree4, node=126, levels_back=0)
cladenrDNA.bs <- data.frame(as_tibble(cladenrDNA))
cladenrDNA.bs$bsvalue <- as.numeric(cladenrDNA.bs$label)
cladenrDNA.bs.node <- subset(cladeplastome.bs, bsvalue >= 90)$node 


T4 <-ggtree(cladenrDNA) + 
  geom_point2(aes(subset=(node %in% cladenrDNA.bs.node)),color="red",size=2) +
  geom_tiplab(aes(color = group.y), align = T, show.legend = F, hjust = 1) +
  scale_color_manual(values=c("G1"= "#D55E00",
                              "G2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "Outgroup" = "black")) +
  geom_treescale(x=0, y =-0.5, fontsize=6, linesize=2, offset=-1.5) +
  scale_y_continuous(expand=c(0, 1.5)) +
  hexpand(1) +
  scale_x_reverse()  


pdf("FigS11_Azorella_plastome_nrDNA3.pdf", width = 10, height = 10)
#fig4
ggarrange(T3,T4, 
          labels = c("a) Plastome - NZ2_cp1", 
                     "b) nrDNA - NZ2"), 
          hjust = -0.5)
dev.off()


png("FigS11_Azorella_plastome_nrDNA3.png", width = 2300, height = 2300, res =250)
ggarrange(T3,T4, 
          labels = c("a) Plastome - NZ2_cp1", 
                     "b) nrDNA - NZ2"), 
          hjust = -0.5)
dev.off()
