setwd("C:/Users/7827x/OneDrive/Desktop/ThesisData/ThesisWriting/Weixuan_Thesis_Chapter2/03_TargetEnrichmentData/02_123samples/03_hypiper")

library("patchwork")     
library("png")
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
library(ggpubr)

Azorella123 <- read.tree(file = "loci_gt07_n90_Azorella.bs30.astral.outgroup.rerooted.treefile")

#############################################################
#############################################################
#############################################################
#############################################################
#############################################################
#############################################################

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
###################black,  orange,    sky blue,  blulsh green, yellow, blue, vermilion, reddish purple 
scale_colour_manual(values=cbPalette)


tax <- read.table("colortable.csv", sep = ',', header  = T)
colnames(tax) = c("treename", "species","accession", "ploidy", "Group")
tax1 <- as.data.frame(tax)


tree1 = full_join(as_tibble(Azorella123), tax1, by = c('label' = 'treename'))
tree1$label2 <- paste0( "[", tree1$ploidy, "] ", tree1$label)
tree2 = as.treedata(tree1)

ggtree(tree2)+   geom_text(aes(label=label2)) + geom_treescale(x=0, y=45)

p <- ggtree(tree2) + 
  geom_tiplab(aes(label= label), align = T, size =8, offset = 0.2,show.legend=FALSE) + 
  geom_nodelab(size = 3, nudge_x = 0.11, color = 'red')  +
  geom_rootedge(rootedge = 0.1)  +
  geom_strip("A.cyanopetala_Lux9", "A.exigua_Dun2", "NZ2", color = "#0072B2", 
             offset = 2.3, align = TRUE, offset.text = 0.05, 
             barsize = 3, fontsize = 10, parse = TRUE) +

  geom_strip("A.fragosea_CANB797854", "A.fragosea_CANB798456", "Au", color = "#CC79A7", 
             offset = 2.3, align = TRUE, offset.text = 0.05, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  
  geom_strip("A.allanii_Rau10", "A.hookeri_AK209018", "NZ1", color = "#D55E00", 
             offset = 2.3, align = TRUE, offset.text = 0.05, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  
  geom_strip("A.burkartii_NYBG2434", "A.ranunculus_NYBG2447", "SA", color = "#000000", 
             offset = 2.3, align = TRUE, offset.text = 0.05, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  geom_strip("A.polaris_End8", "A.robusta_CHR617278", "Sub", color = "#009E73", 
             offset = 2.3, align = TRUE, offset.text = 0.05, 
             barsize = 3, fontsize = 10, parse = TRUE) + 
  geom_treescale(x=0, y=-1, linesize = 3, fontsize=6, offset= -1)
p


#############################################################
#############################################################
#############################################################

k <- read.table("out.node.key", header=F, row.names=NULL, sep=" ", na.strings = c("", "NA"), stringsAsFactors=FALSE)

colnames(k) <- c("NodeP", "Tips")
NodeP <- rep(0, length(k$NodeP))
k2 <- as.data.frame(NodeP)
k2$NodeT <- rep(0, length(k2$NodeP))

for (x in 1:length(k[,1])) {
  k2$NodeP[x] <-  k[x,1]
  tips <- gsub("\\)|\\(", ",", k[x,2])
  tips_vec <- strsplit(tips, ",")
  tips_vec2 <- grep("A.*",tips_vec[[1]], value=T)
  k2$NodeT[x] <- mrca.phylo(Azorella123, tips_vec2)
}


phypart_pie <- read.table(file = "phyparts_pies.csv", sep = ",", header = T)
new_pie <- cbind(k2, phypart_pie)[c(2,4:7)]
names(new_pie)[names(new_pie) == 'NodeT'] <- 'node'

pies <- nodepie(new_pie, cols = 2:5, alpha = 0.8,
                color = c(adj_concord = "#56B4E9", 
                         adj_most_conflict = "#009E73", 
                         other_conflict = "#CC79A7", 
                         the_rest = "grey"))


p1 <- p  + geom_inset(pies, width = 0.06, height = 0.06, hjust = 0.02,  vjust = 0.06) 
p1

#############################################################
#############################################################
#############################################################

sample.data= as.matrix(read.csv(file = "gene_lengths123.txt",header=T,row.names=1,sep="\t"))
sample.len = sample.data[2:nrow(sample.data),]
reference.len = as.numeric(sample.data[1,])

#Calculate the percentage length recovered relative to the reference.
percent.len=sweep(sample.len,2,as.numeric(reference.len),'/')
percent.len = ifelse(percent.len>1,1,percent.len)

percent.long = melt(percent.len)

#############################################################
#############################################################
#############################################################
#snp.data = as.matrix(read.csv(file = "0_Table_SNPs.csv",header=T,row.names=1,sep=","))
#write.csv(t(snp.data), file = "SNP_boxplot.csv", quote = F)

snp.data.renamed = as.matrix(read.csv(file = "SNP_boxplot.csv",header=T,row.names=1,sep=","))
snp.data.renamed_melt <- melt(snp.data.renamed)
snp.data.renamed_melt$Var1
#############################################################
#############################################################
#############################################################

p2 <- p1 + geom_fruit(data=percent.long, geom=geom_tile,
                     mapping= aes(x=Var2, y=Var1, fill = value),
                     offset = 0.31, pwidth = 0.6, size = 3,
                     axis.params=list( axis = "x", text.angle = 90, size = 0.01)) +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_gradient(high = "#132B43", low = "#FFFFFF")

p2

cols <- c("?x" = "grey", "2x" = "#CC79A7", "4x" = "skyblue", "6x" = "seagreen3", "10x" = "yellow3")

p3 <- p2  + new_scale_fill() +  
  geom_fruit( data=snp.data.renamed_melt,
              geom=geom_boxplot,
              mapping = aes(y=Var1, x = value, color = 'black', fill = ploidy),
              size=0.2, offset = 0.06, pwidth = 0.45,
              outlier.size=0.5,   outlier.stroke=0.2,
              outlier.shape=21, 
              axis.params=list( axis = "x",  text.size  = 6,
                                hjust = 0.5, vjust = 0.9, 
                                nbreak = 6),
              grid.params=list(size = 0.5)) +
  scale_fill_manual(values=cols) +
  theme(legend.position = c(0.98, 0.5),
        legend.key.size = unit(1, 'cm'),
        legend.title=element_text(size=14),
        legend.text=element_text(size= 20)) 
  
p3
#############################################################
#############################################################
#############################################################
library(scales)
library(ggrepel)

slices <- c(90, 90, 90, 90)
lbls <- c("concordant", "discordant with\n a main alternative", "discordant with\n all remaining alternatives", "uninformative")

pie.data <- data.frame(group = lbls, value = slices)

cols <- c("#56B4E9",  "#009E73",  "#CC79A7",  "grey")

pie.data <- pie.data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(pie.data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )


png("pie_char.png", width = 3000, height = 3000, res =200)

ggplot(pie.data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=0.1, color="white", alpha = 0.8, show.legend=FALSE) +
  coord_polar("y", start=0) +
  scale_fill_manual(values=cols) +
  theme_void() 

dev.off()
#############################################################
#############################################################

#png("Azorella_targeted123.png", width = 4000, height = 6000, res =250)
#gheatmap(p, percent.len, offset=5.5, width=5, font.size=1, 
#         colnames_angle=-45, hjust=0) +
#  scale_fill_gradient(high = "#132B43", low = "#FFFFFF")
#dev.off()




png("Azorella_targeted_heatmap_SNP.png", width = 12000, height = 15000, res =400)
p3
dev.off()

pdf("Fig3_Azorella_targeted_heatmap_SNP.pdf", width = 30, height = 40)
p3
dev.off()


#############################################################
#############################################################
#############################################################
my_image <- readPNG("pie_char.png", native = TRUE)

ggp_image <- p3 +                  # Combine plot & image
  inset_element(p = my_image,
                left = 0.01,
                bottom = 0.75,
                right = 0.21,
                top = 0.85)

pdf("Fig3A_Azorella_targeted_heatmap_SNP.pdf", width = 30, height = 40)
ggp_image
dev.off()

png("Fig3A_Azorella_targeted_heatmap_SNP.png", width = 12000, height = 15000, res =400)
ggp_image
dev.off()
#############################################################
#############################################################
#############################################################

Azorellaconcat <- read.tree(file = "concord.cf.tree")
Azorellaconcat2 <- root(Azorellaconcat, outgroup = "A.lycopodioides_NYBG2433", edgelabel = TRUE)

p4 <- ggtree(Azorellaconcat2) + 
  geom_tiplab(align = T, size =8, offset = 0.01) + 
  geom_nodelab(size = 3, nudge_x = 0.001, color = 'red') +
  xlim(c(0,0.2)) +
  geom_strip("A.pallida_Cra4", "A.pallida_CHR501035", "NZ2", color = "#0072B2", 
             offset = 0.1, align = TRUE, offset.text = 0.005, 
             barsize = 3, fontsize = 10, parse = TRUE) +

  
  geom_strip("A.fragosea_CANB798456", "A.fragosea_CANB797887", "Au", color = "#CC79A7", 
             offset = 0.1, align = TRUE, offset.text = 0.005, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  
  geom_strip("A.allanii_HikA9",  "A.nitens_CHR252408", "NZ1", color = "#D55E00", 
             offset = 0.1, align = TRUE, offset.text = 0.005, 
             barsize = 3, fontsize = 10, parse = TRUE) +

  
  geom_strip("A.burkartii_NYBG714" , "A.ranunculus_NYBG2447", "SA", color = "#000000", 
             offset = 0.1, align = TRUE, offset.text = 0.005, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  geom_strip("A.polaris_End8", "A.lyallii_CHR542359", "Sub", color = "#009E73", 
             offset = 0.1, align = TRUE, offset.text = 0.005, 
             barsize = 3, fontsize = 10, parse = TRUE) +
  geom_treescale(x=0, y=-1, linesize = 3, fontsize=6, offset= -1) +
  #ggtitle("123 Sample Concatenated Tree") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = -1)) 

p4


png("FigS4_Azorella_targeted_concat.png", width = 5000, height = 10000, res =300)

p4
dev.off()


pdf("FigS4_Azorella_targeted_concat.pdf", width = 20, height = 40)

p4
dev.off()

#######################################################################
########################################################################
png("Azorella_targeted_heatmap_SNP_concat2.png", width = 16000, height = 20000, res =400)
ggarrange(p3, p4,widths = c(1.2, 0.5),
          labels = c("a)", "b)"), 
          font.label = list(size = 30),
          ncol = 2, nrow = 1)
dev.off()

pdf("Azorella_targeted_heatmap_SNP_concat.pdf", width=23, height=25)
ggarrange(p3, p4,widths = c(1.5, 0.5),
          labels = c("a)", "b)"), 
          font.label = list(size = 30),
          ncol = 2, nrow = 1)
dev.off()
