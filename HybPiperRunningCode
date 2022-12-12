library(scales)
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
library(tanggle)
library(phangorn)

#Setup your working directory here

setwd("C:/Users/7827x/OneDrive/Desktop/Example_LPP_Pie")

#read your ASTRAL tree loci_gt07_n4.astral.rerooted.treefile

Azorella.illumina <- read.tree(file = "loci_gt07_n4.astral.rerooted.treefile")


##########################################################################################

# This section is to use the phyparts output as input to plot on the tree

# Phyparts has two outputs 1) phyparts_pies.csv that calculated for the proportion of each pie chart
# 2) phyparts_dist.csv that calcluate the number of concordant gene vs all the remaining sample. 


# The big challenging here is that the nodes of phyparts output (please have a look at the output file), 
# will not be the same of the nodes in your tree. 

# To assign the pie chart and values to the correct node of the tree, we have to use the function "mrca.phylo" to search for
# the nodes that contain the same tips between phyparts outcome and ASTRAL tree, and merge them together,
# so that, your pie charts will be on the correct position of your tree.


# k below is the nodes of phyparts

k <- read.table("Az_illumina_n4_bs30.node.key", header=F, row.names=NULL, sep=" ", na.strings = c("", "NA"), stringsAsFactors=FALSE)

colnames(k) <- c("NodeP", "Tips")
NodeP <- rep(0, length(k$NodeP))
k2 <- as.data.frame(NodeP)
k2$NodeT <- rep(0, length(k2$NodeP))

# k2 below is the nodes of ASTRAL tree

for (x in 1:length(k[,1])) {
  k2$NodeP[x] <-  k[x,1]
  tips <- gsub("\\)|\\(", ",", k[x,2])
  tips_vec <- strsplit(tips, ",")
  tips_vec2 <- grep("A.*",tips_vec[[1]], value=T)
  k2$NodeT[x] <- mrca.phylo(Azorella.illumina, tips_vec2)
}



##########################################################################################

# Read outputs from phyparts and converted the nodes number to ASTRAL tree corespned number


phypart_pie <- read.table(file = "phyparts_pies.csv", sep = ",", header = T)
new_pie <- cbind(k2, phypart_pie)[c(2,4:7)]
names(new_pie)[names(new_pie) == 'NodeT'] <- 'node'

pies.illumina <- nodepie(new_pie, cols = 2:5, alpha = 0.8,
                         color = c(adj_concord = "#56B4E9", 
                                   adj_most_conflict = "#009E73", 
                                   other_conflict = "#CC79A7", 
                                   the_rest = "grey"))




phypart_pie_node <- read.table(file = "phyparts_dist.csv", sep = ",", header = T)
phypart_pie_node$plot <- str_c(phypart_pie_node$concord, "/",phypart_pie_node$genes.concord)
new_pie_node <- cbind(k2, phypart_pie_node)[c(2,6)]
names(new_pie_node)[names(new_pie_node) == 'NodeT'] <- 'node'




##########################################################################################

# Here, to make the tips colorful, I manunlly wrote a table based on their genetic groups called colortable.csv
# I then added the phyparts output to the phylogenetic tree by first converting the tree into a table "as_tibble"
# then joined the data onto the table "full_join", then converted it back to a tree "as.treedata"

tax <- read.table("colortable.csv", sep = ',', header  = T)
colnames(tax) = c("treename", "species","accession", "ploidy", "Group")
tax1 <- as.data.frame(tax)

test <- as_tibble(Azorella.illumina)
tree1 = full_join(as_tibble(Azorella.illumina), new_pie_node, by = c('node' = 'node'))
tree2 = full_join(tree1, tax1, by = c(label = 'treename'))

tree2$label2 <- paste0( "[", tree2$ploidy, "] ", tree2$label)

tree3 = as.treedata(tree2)

#####################

p.iilumina <- ggtree(tree3) + 
  geom_text(aes(label=plot), size = 3, vjust =  -0.2, hjust= 1.5) + 
  geom_tiplab(aes(label = label, color = Group), align = T, size =5, offset = 0.1, show.legend=FALSE) + 
  geom_nodelab(size = 3, nudge_x = 0.3, vjust= 0.5, color = 'red')  +
  geom_rootedge(rootedge = 0.1) +
  #  ggtitle("ASTRAL Tree of 22 Taxa") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, vjust = -1)) +
  scale_color_manual(values=c("G1"= "#D55E00",
                              "G2"= "#0072B2",
                              "SUB"= "#009E73",
                              "AUZ" = "#CC79A7",
                              "Outgroup" = "black")) +
  xlim(0, 8) +
  geom_treescale(x=0, y =-0.05, fontsize=3, linesize=2, offset=-0.5) 

p2.iilumina <- p.iilumina   + geom_inset(pies.illumina, width = 0.15, height = 0.15, hjust = 0.02,  vjust = 0.06) 

p2.iilumina

####################################################

#Addd the heatmap from Hybpiper2 output using file seq_length.tsv, make sure your sample names are identical to your tree tip names

sample.data= as.matrix(read.csv(file = "seq_lengths.tsv",header=T,row.names=1,sep="\t"))
sample.len = sample.data[2:nrow(sample.data),]
reference.len = as.numeric(sample.data[1,])

percent.len =sweep(sample.len,2,as.numeric(reference.len),'/')
percent.len = ifelse(percent.len>1, 1, percent.len)
percent.long = melt(percent.len)


p3.iilumina <- p2.iilumina + geom_fruit(data=percent.long, geom=geom_tile,
                      mapping= aes(x=Var2, y=Var1, fill = value),
                      offset = 0.4, pwidth = 0.7, size = 4,
                      axis.params=list( axis = "x", text.angle = 90, size = 0.01)) +
  scale_fill_gradient(high = "#132B43", low = "#FFFFFF") +
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = c(0.97, 0.5),
        legend.key.size = unit(0.7, 'cm'),
        legend.title=element_text(size=10),
        legend.text=element_text(size= 10)) 

p3.iilumina
