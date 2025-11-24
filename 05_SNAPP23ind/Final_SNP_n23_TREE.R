setwd("C:/Users/7827x/OneDrive/Desktop/Azorella_NewBatch/01_Chapter2/01_Angiosperms353/06_n23_timetree")

library(ggtree)
library(treeio)
library(tidytree)
library(dplyr)
library(ggplot2)
library(aplot)
library(ggtreeExtra)
library(ggnewscale)
library(reshape2)
library(phangorn)
library(stringr)
library(ggpubr)
library(phangorn)
library("patchwork")     
library("png")

#remotes::install_github("KlausVigo/ggnetworx")
#install.packages("deeptime")
library(ggnetworx)
library(deeptime)
library(phytools)

#######################

snptree <- read.beast.newick(file = "Azorella_n23.snp.filtered.rename.bialleic.thin20.min4_reroot.tree")

ggtree(snptree) +
  geom_tiplab(aes(label = label.x)) +
  geom_nodelab(aes(label = label.y), size = 4, nudge_x = 0.28,  color = 'red')  


tax <- read.table("colortable2.csv", sep = ',', header  = T)
tree1 = full_join(as_tibble(snptree), tax, by = c(label.x = 'TreeName'))
colnames(tree1)[4] <- "label"
colnames(tree1)[5] <- "bbvalue"

tree2 = as.treedata(tree1)


ggtree(tree2)

svdtree <- ggtree(tree2) + 
  geom_tiplab(aes(color = Group),  size =6, offset = 0.2, show.legend=FALSE) + 
  geom_nodelab(aes(label = bbvalue), size = 4, nudge_x = 0.28,  color = 'black')  +
  geom_rootedge(rootedge = 0.2) +
  
  #  geom_text(aes(label=round(posterior,2)), size = 3, hjust= -0.5) +
  
  #  ggtitle("SVDQuartets SNP Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, vjust = -1)) +
  geom_treescale(x=0.1, fontsize=6, linesize=1, offset = -0.5)  +
  scale_y_continuous(expand=c(0, 1)) +
  scale_color_manual(values=c("G1"= "#D55E00",
                              "G2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "Outgroup" = "black")) +
  xlim(-0.2,15) 

svdtree


###################################################################################
###################################################################################

snappptree21 <- read.beast("snapp.tre")
ggtree(snappptree21) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab() 


tax <- read.table("colortable2.csv", sep = ',', header  = T)
colnames(tax) = c("treename", "species","accession", "ploidy", "Group")
tax1 <- as.data.frame(tax)

tree21.snapp = full_join(as_tibble(snappptree21), tax1, by = c(label = 'treename'))
tree21.snapp2 = as.treedata(tree21.snapp)



ggtree(tree21.snapp2) + geom_text(aes(label=node))

heighthpd <- as_tibble(snappptree21)
heighthpd.sub <- heighthpd[c(24,25,26,43,28,30,27,33, 34, 35),]
newtable <- data.frame(round(t(data.frame(heighthpd.sub$height_0.95_HPD)),2))
newtable$hpd <- round(heighthpd.sub$height,2)
newtable$node <- heighthpd.sub$node
newtable$hpdplot <- paste0(newtable$hpd, " [", newtable$X2, ", ", newtable$X1, "]" )

nodeplot <- data.frame(newtable$node,newtable$hpdplot)

tree21.snapp3 = full_join(as_tibble(tree21.snapp2), nodeplot, by = c(node = 'newtable.node'))
tree21.snapp4 = as.treedata(tree21.snapp3)



p2 <- ggtree(tree21.snapp4) + 
  geom_text(aes(label=newtable.hpdplot), size = 3, hjust= 1.2, vjust=-1) +
  geom_text(aes(label=round(posterior,2)), size = 3, hjust= -0.5,  color = "blue") +
  geom_tiplab(aes(color = Group),  size =6, offset = 1, align = T, show.legend=FALSE) +
  scale_color_manual(values=c("G1"= "#D55E00",
                              "G2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "Outgroup" = "black")) +
  theme_tree2(axis.text.x = element_text(size =20),
              axis.line.x = element_line(colour = 'black', size = 2),
              axis.ticks.x = element_line(colour = "black", size = 2),
              plot.title = element_text(hjust = 0.5, size = 18, vjust = -1))
p2

snappptree21.2 <- revts(p2) +  
  geom_rootedge(rootedge = 1) +
  scale_x_continuous(limits = c(-37, 12), breaks=seq(-35, 0, 5), labels=abs) +
  scale_y_continuous(expand=c(0, 2)) +
  geom_range('height_0.95_HPD', color='red',  size=3, alpha=.3)  +  
  xlab("Time Ma")  +
  theme(axis.title=element_text(size=10,face="bold", hjust=0.7)) 
#+  ggtitle("Baysiane SNAPP Tree of 21 Taxa") 

snappptree21.2


###################################################################################
###################################################################################

phypart_pie <- read.table(file = "range_probabilities.txt", sep = "\t", header = T)
phypart_pie$anslabel <- colnames(phypart_pie)[max.col(phypart_pie,ties.method="first")]
phypart_pie$node <- as.numeric(rownames(phypart_pie))

##########cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
###################black,  orange,    sky blue,  blulsh green, yellow, blue, vermilion, reddish purple 

fill.biogeo <- c(A = "#56B4E9", B = "#009E73", C = "#CC79A7", D = "#D55E00", E = "#0072B2", 
                 "F" = "#E69F00", "DE" = "#B2182B", "AB" = "#0000FF","EF" = "#000000")
pies <- nodepie(phypart_pie[c(24:45),], cols = 1:22, alpha = 0.8, 
                color = c(A = "#56B4E9", 
                          B = "#009E73", 
                          C = "#CC79A7", 
                          D = "#D55E00",
                          E = "#0072B2",
                          "F" = "#E69F00",
                          "DE" = "#B2182B",
                          "AB" = "#0000FF",
                          "EF" = "#000000"))

snappptree <- read.beast.newick("../08_n23_BioGeo/21snapp.tre.nex") 
tree.ans = full_join(as_tibble(snappptree), phypart_pie[23:24], by = c('node' = 'node'))

tree.ans$treename <- paste0("[", tree.ans$anslabel,"] ", tree.ans$label )
tree.ans.snapp = as.treedata(tree.ans)


p <- ggtree(tree.ans.snapp) + 
  geom_text2(aes(subset=!isTip, label=anslabel, col = as.factor(anslabel)), hjust=1.7, vjust =-1.7, show.legend=FALSE) + 
  geom_tiplab(aes(color = as.factor(anslabel), label = treename),  size =6, offset = 1, align = T, show.legend=FALSE) +
  scale_color_manual(values = fill.biogeo) +
  theme_tree2(axis.text.x = element_text(size =20),
              axis.line.x = element_line(colour = 'black', size = 2),
              axis.ticks.x = element_line(colour = "black", size = 2),
              plot.title = element_text(hjust = 0.5, size = 18, vjust = -1)) 

p
snappptree.biogeo <- revts(p) +  
  geom_rootedge(rootedge = 1) +
  scale_x_continuous(limits = c(-30, 13), breaks=seq(-30, 0, 5), labels=abs) +
  scale_y_continuous(expand=c(0, 2)) +
  xlab("Time Ma")  +
  
  geom_strip("Azpallida_Cra4", "Azhaastii_Eri1", "NZ2", color = "#0072B2", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) +
  geom_strip("Azschizeilema_Hoo", "Azschizeilema_Hoo", "Sub", color = "#009E73", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) + 
  geom_strip("Azhaastii_Bri6", "Azexigua_Wet1", "NZ2", color = "#0072B2", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) +
  geom_strip("Azfragosea_CANB797854", "Azfragosea_CANB797854", "Au", color = "#CC79A7", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) +
  
  
  geom_strip("Azhookeri_Tai", "Azallanii_HikB1", "NZ1", color = "#D55E00", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) +
  
  geom_strip("Azranunculus_NYBG2447", "Azburkartii_NYBG714", "SA", color = "black", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) +
  
  geom_strip("Azpolaris_Camp2", "Azrobusta_CHR617278", "SubM", color = "#009E73", 
             offset = 12.2, align = TRUE, offset.text = 0.5, 
             barsize = 4, fontsize = 6, parse = TRUE) + 
  
  annotate("segment", x = 12.2, xend = 12.2, y = 12.7, yend = 13.2, size =4, colour = "#CC79A7") +
  annotate("segment", x = 12.2, xend = 12.2, y = 16.7, yend = 17.2, size =4, colour = "#009E73") +
  
  theme(axis.title=element_text(size=10,face="bold", hjust=0.7)) 
#+  ggtitle("Ancestral Range Reconstruction with BioGeoBear for 21 Taxa") 

snappptree.biogeo
snappptree.biogeo.pie <- snappptree.biogeo  + geom_inset(pies, width = 0.08, height = 0.08, hjust = 0.02,  vjust = 0.06) 

snappptree.biogeo.pie.legend <- snappptree.biogeo.pie + 
  annotate("rect", xmin = -25.5, xmax = -14, ymin = 11, ymax = 21, alpha = 0.1) +
  annotate("text", x = -25, y = 20.5, size =5, label = "A = South America (SA)", color = "#56B4E9", hjust = 0) + 
  annotate("text", x = -25, y = 19.5, size =5, label = "B = Subantarctic Islands (Sub), NZ", color = "#009E73", hjust = 0) +
  annotate("text", x = -25, y = 18.5, size =5, label = "C = Australia", color = "#CC79A7", hjust = 0) +
  annotate("text", x = -25, y = 17.5, size =5, label = "D = Nouth Island (N), NZ", color = "#D55E00", hjust = 0) +
  annotate("text", x = -25, y = 16.5, size =5, label = "E = Sorth Island (S), NZ", color = "#0072B2", hjust = 0) +
  annotate("text", x = -25, y = 15.5, size =5, label = "F = Stewart Island", color = "#E69F00", hjust = 0) +
  annotate("text", x = -25, y = 14.5, size =5, label = "AB = SA & Sub", color = "#0000FF", hjust = 0) +
  annotate("text", x = -25, y = 13.5, size =5, label = "DE = NZ N/S ", color = "#B2182B", hjust = 0) +
  annotate("text", x = -25, y = 12.5, size =5, label = "EF = NZ S & Stewart Island", color = "#000000", hjust = 0) +
  annotate("text", x = -25, y = 11.5, size =5, label = "All less probable ranges", color = "grey40", hjust = 0) 


snappptree.biogeo.pie.legend +
  geom_text(aes(label=round(posterior,2)), size = 3, hjust= -0.5,  color = "blue") 



###################################################################################
###################################################################################
###################################################################################
g <- ggarrange(svdtree, NULL,  snappptree21.2, snappptree.biogeo.pie.legend,
               labels = c(#"a) MCMC_BiMarkers Network",
                 "a) SVDquartets ML tree",
                 "b) Cloudogram of SNAPP Trees", 
                 "c) Consensus SNAPP Tree", 
                 "d) Ancestral Range Reconstruction"), 
               hjust = -0.05,
               font.label = list(size = 20),
               ncol = 2, nrow = 2)
g

###################################################################################
###################################################################################

library(RColorBrewer)
library(R.utils)
library(ggplotify)

setwd("C:/Users/7827x/OneDrive/Desktop/Azorella_NewBatch/01_Chapter2/01_Angiosperms353/06_n23_timetree/Treemixoutput")
source("plotting_funcs.R") # here you need to add the path

#prefix="TreeMix.2"

#pdf("test.pdf", height = 10, width = 10)
#par(mfrow=c(2,3))
#for(edge in 1:5){
#  plot_tree(cex=0.8,paste0(prefix,".",edge))
#  title(paste(edge,"edges"))
#}
#dev.off()

plot_tree(cex=1,  arrow = 0.1, scale = F, plu = 0.05, lwd = 2,  mbar  = F, "TreeMix.2.1")

test <- as.ggplot(~plot_tree(cex=1.4, font =1, arrow = 0.1, scale = F, plu = 0.05, lwd = 2,  "TreeMix.2.1"))

g <- ggarrange(test, 
               NULL,  snappptree21.2, snappptree.biogeo.pie.legend,
               labels = c(#"a) MCMC_BiMarkers Network",
                 "a) TreeMix (m=2)",
                 "b) Cloudogram of SNAPP Trees", 
                 "c) Consensus SNAPP Tree", 
                 "d) Ancestral Range Reconstruction"), 
               hjust = -0.05,
               font.label = list(size = 20),
               ncol = 2, nrow = 2)


setwd("C:/Users/7827x/OneDrive/Desktop/Azorella_NewBatch/01_Chapter2/01_Angiosperms353/06_n23_timetree")

###################################################################################
###################################################################################
my_image <- readPNG("23SNAPP_tree_noname_grid.png", native = TRUE)

Fig8 <- g +                
  annotate("text", x = 0.805, y = 0.96, size = 6, label = "[4x] Azschizeilema_Hoo", color = "#009E73",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.94, size = 6, label = "[4x] Azhaastii_Bri6", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.92, size = 6, label = "[6x] Azpallida_Bro4", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.903, size = 6, label = "[6x] Azpallida_Cra4", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.885, size = 6, label = "[4x] Azhaastii_Pat10A", color = "#0072B2",   hjust = 0) +
  annotate("text", x = 0.805, y = 0.865, size = 6, label = "[4x] Azcyanopetala_Tak5", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.846, size = 6, label = "[4x] Azcyanopetala_Bro3", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.827, size = 6, label = "[4x] Azhaastii_Eri1", color = "#0072B2", hjust = 0) +
  annotate("text", x = 0.805, y = 0.806, size = 6, label = "[4x] Azhydrocotyloides_Dob6", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.787, size = 6, label = "[4x] Azexigua_Wet1", color = "#0072B2",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.77, size = 6, label = "[?x] Azfragosea_CANB797854", color = "#CC79A7",  hjust = 0) +
  
  annotate("text", x = 0.805, y = 0.75, size = 6, label = "[6x] Azhookeri_Tai", color = "#D55E00",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.73, size = 6, label = "[10x] Azcolensoi_BellA10", color = "#D55E00",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.71, size = 6, label = "[6x] Azcockaynei_Man10", color = "#D55E00", hjust = 0) +
  annotate("text", x = 0.805, y = 0.692, size = 6, label = "[6x] Aznitens_Wai3", color = "#D55E00", hjust = 0) +
  annotate("text", x = 0.805, y = 0.672, size = 6, label = "[4x] Azroughii_Sta1", color = "#D55E00",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.655, size = 6, label = "[4x] Azallanii_HikB1", color = "#D55E00",   hjust = 0) +
  
  annotate("text", x = 0.805, y = 0.635, size = 6, label = "[?x] Azburkartii_NYBG714", color = "black",  hjust = 0) +
  annotate("text", x = 0.805, y = 0.615, size = 6, label = "[2x] Azranunculus_NYBG2447", color = "black",  hjust = 0) +
  
  annotate("text", x = 0.805, y = 0.595, size = 6, label = "[6x] Azpolaris_Camp2", color = "#009E73", hjust = 0) +
  annotate("text", x = 0.805, y = 0.576, size = 6, label = "[?x] Azlyallii_CHR542359", color = "#009E73", hjust = 0) +
  annotate("text", x = 0.805, y = 0.56, size = 6, label = "[?x] Azrobusta_CHR617278", color = "#009E73", hjust = 0) +
  
  annotate("text", x = 0.805, y = 0.54, size = 6, label = "[2x] Azlycopodioides_NYBG2433", color = "black", hjust = 0) +
  #  annotate("text", x = 0.79, y = 0.54, size = 6, label = "[10x] Azcolensoi_BellA10", color = "#D55E00", hjust = 0) +
  
  annotate("segment", x = 0.935, xend = 0.935, y = 0.955, yend = 0.965, size =4, colour = "#009E73") +
  annotate("segment", x = 0.935, xend = 0.935, y = 0.787, yend = 0.94, size =4, colour = "#0072B2") +
  annotate("segment", x = 0.935, xend = 0.935, y = 0.765, yend = 0.775, size =4, colour = "#CC79A7") +
  annotate("segment", x = 0.935, xend = 0.935, y = 0.65, yend = 0.75, size =4, colour = "#D55E00") +
  annotate("segment", x = 0.935, xend = 0.935, y = 0.61, yend = 0.64, size =4, colour = "black") +
  annotate("segment", x = 0.935, xend = 0.935, y = 0.56, yend = 0.59, size =4, colour = "#009E73") +
  
  annotate("text", x = 0.94, y = 0.96, size = 6, label = "Sub", color = "#009E73", hjust = 0) +
  annotate("text", x = 0.94, y = 0.85, size = 6, label = "NZ2", color = "#0072B2", hjust = 0) +
  annotate("text", x = 0.94, y = 0.77, size = 6, label = "Au", color = "#CC79A7", hjust = 0) +
  annotate("text", x = 0.94, y = 0.7, size = 6, label = "NZ1", color = "#D55E00", hjust = 0) +
  annotate("text", x = 0.94, y = 0.625, size = 6, label = "SA", color = "black", hjust = 0) +
  annotate("text", x = 0.94, y = 0.575, size = 6, label = "SubM", color = "#009E73", hjust = 0) +
  
  
  inset_element(p = my_image,
                on_top = FALSE,
                left = 0.5,
                bottom = 0.52,
                right = 0.995,
                top = 0.97)

Fig8

###################################################################################
###################################################################################
png("Final_Fig7redo_Azorella22_SNPTree.png", width = 8000, height = 5500, res =300)
Fig8
dev.off()

pdf("Final_Fig7redo_Azorella22_SNPTree.pdf", width = 28, height = 20)
Fig8
dev.off()

