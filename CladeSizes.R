library(castor)
library(geiger)
library(dispRity)
library(TreeTools)

#Load a tree from a string or file in Newick (parenthetic) format
tree<-read_tree(file="phrynosomatid.tre")

#Clade sizes
CladeSizes(tree, nodes = c(114, 117, 9))

#Tips contained within splits
splits <- as.Splits(tree)
summary(splits)
TipsInSplits(splits)

#extract.clade
extract_clade<-extract.clade(tree, node=117, root.edge = 0, collapse.singles = TRUE,
                             interactive = FALSE)
extract_clade
write.nexus(extract_clade,file="subtree117")


#
CladisticInfo(tree)
ct1 <- as.ClusterTable(tree)
summary(ct1)
as.matrix(ct1)

#Convert object to Splits
tree <- BalancedTree(LETTERS[1:5])
splits <- as.Splits(tree)
plot(tree)
summary(splits)
LabelSplits(tree, as.character(splits), frame = "none", pos = 3L)
LabelSplits(tree, TipsInSplits(splits), unit = " tips", frame = "none",
            pos = 1L)
splits1 <- as.Splits(BalancedTree(7))
splits2 <- as.Splits(PectinateTree(7))
match(splits1, splits2)

#Count descendants for each node in a tree
NDescendants(tree)
nodelabels(NDescendants(tree))

#Subsplit
splits <- as.Splits(PectinateTree(letters[1:9]))
splits
efgh <- Subsplit(splits, tips = letters[5:8], keepAll = TRUE)
summary(efgh)
TrivialSplits(efgh)
summary(Subsplit(splits, tips = letters[5:8], keepAll = FALSE))


TipLabels(tree)

splits <- as.Splits(tree)
TipsInSplits(splits)
LabelSplits(tree, as.character(splits))
as.TreeNumber(tree)

splitTree(tree, split)
