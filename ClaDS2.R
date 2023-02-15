library(ape)
library(castor)
library(geiger)
library(dispRity)
library(phytools)
library(caper)
library(TreeTools)

#Load a tree from a string or file in Newick (parenthetic) format
tree<-read_tree(file="Archaea_species.nwk")
write.tree(tree,file="Archaea_species.tre")
#Extracting the age of nodes in a tree
node_age=branching.times(tree)
write.csv(node_age,file="node_age")

#List the node index
node_index=c(Ntip(tree)+1:Nnode(tree))
write.csv(node_index,file="node_index")

# Create a empty matrix
mat = matrix(, nrow = Nnode(tree), ncol = 8)

# Put the species richness of each node in the matrix
for (x in (1:Nnode(tree))) {
  mat[x,1] <- node_age[x]#the age of nodes
  mat[x,2] <- length(tips(tree,node_index[x]))#species richness
  mat[x,3] <- log10(length(tips(tree,node_index[x])))
  mat[x,4] <- node_index[x]#node index of pick up clades
  mat[x,5] <- bd.ms(phy=extract.clade(tree, 
                                      node=node_index[x]),
                    missing=0, epsilon=0)
  mat[x,6] <- bd.ms(phy=extract.clade(tree, 
                                      node=node_index[x]),
                    missing=0, epsilon=0.5)
  mat[x,7] <- bd.ms(phy=extract.clade(tree, 
                                      node=node_index[x]),
                    missing=0, epsilon=0.9)
  combindedtips <- list(tips(tree,node_index[x]))
  mat[x,8] <- paste0(combindedtips,sep=",")
}
#print(mat)
write.csv(mat,file="mat")


#Randomly pick up data
for (x in (1:10)){
  df_sample = mat[sample(nrow(mat), 50), ]
  filename1 <- paste0("Randomly pick up data_", x)
  
  #name the columns
  colnames(df_sample) <- c('Clage_age','SR',
                           'Log10SR',
                           'node_index','div_rate_e0',
                           'div_rate_e5','div_rate_e9',
                           'species')
  write.csv(df_sample,file=filename1)
}

#combine the results
Rep1 <- read.csv(file="Randomly pick up data_1",
                 stringsAsFactors = TRUE)
Rep2 <- read.csv(file="Randomly pick up data_2",
                 stringsAsFactors = TRUE)
Rep3 <- read.csv(file="Randomly pick up data_3",
                 stringsAsFactors = TRUE)
Rep4 <- read.csv(file="Randomly pick up data_4",
                 stringsAsFactors = TRUE)
Rep5 <- read.csv(file="Randomly pick up data_5",
                 stringsAsFactors = TRUE)
Rep6 <- read.csv(file="Randomly pick up data_6",
                 stringsAsFactors = TRUE)
Rep7 <- read.csv(file="Randomly pick up data_7",
                 stringsAsFactors = TRUE)
Rep8 <- read.csv(file="Randomly pick up data_8",
                 stringsAsFactors = TRUE)
Rep9 <- read.csv(file="Randomly pick up data_9",
                 stringsAsFactors = TRUE)
Rep10 <- read.csv(file="Randomly pick up data_10",
                  stringsAsFactors = TRUE)

result <- rbind(Rep1,Rep2,Rep3,Rep4,Rep5,
                Rep6,Rep7,Rep8,Rep9,Rep10)
write.csv(result,file="result")


#read the sampled data

rep <- read.csv(file="result")

#get the picked nodes
nodes = rep$node_index-Ntip(tree)
# nodes = subset(nodes, rep$SR!=2)
#get the tips that included in 50 picked nodes
subtrees = get_subtrees_at_nodes(tree, nodes)$subtrees

for (x in (1:500)){
  filename1 <- paste0(x,".tre")
  write.tree(subtrees[[x]],filename1)
  # if (Ntip(subtrees[[x]])==2){
  #   file.remove(file=filename1)
  # }
}

rate<-data.frame()
for (x in (14:19)){
  if (Ntip(subtrees[[x]])==2){
    rate <- rbind (rate, 0)
  }else if (Ntip(subtrees[[x]])>2) {
    filename1 <- paste0(x)
    load(file=filename1)
    rate <- rbind (rate, CladsOutput$lambdai_map[1])
  }
}
write.csv(rate,file="rate.csv")

rate<-CladsOutput$lambdai_map
write.csv(rate,file="tree.csv")
edge<-tree$edge
write.csv(edge,file="edge.csv")
CladsOutput$lambda0_map
CladsOutput$lambdatip_map

load(file="6")
rate<-CladsOutput$lambdai_map
write.csv(rate,file="6.csv")
edge<-tree1$edge
write.csv(edge,file="6edge.csv")

load(file="12")
rate<-CladsOutput$lambdai_map
write.csv(rate,file="12.csv")
edge<-tree1$edge
write.csv(edge,file="12edge.csv")

CladsOutput$lambdai_map[1]


load(file="13")
rate<-CladsOutput$lambdai_map
write.csv(rate,file="13.csv")
edge<-tree1$edge
write.csv(edge,file="13edge.csv")