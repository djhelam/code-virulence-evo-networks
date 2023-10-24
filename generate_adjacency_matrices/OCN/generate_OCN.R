#generate OCNs used in the paper
library(OCNet) #package by Carraro et al. (2020, Ecology and Evolution)
library(igraph)
source("../functions.R")
set.seed(100) #set seed for random number genertaor
REPLICATES=1000 #number of realisations
i=1
  system(paste("mkdir ",i,sep=""))
  system(paste("mkdir ",i,"/matrices",sep=""))
  for(r in 0:(REPLICATES-1))
{
    #unaggregated OCN on a 10x10 grid with fixed outlet position
        ocn=create_OCN(10,10, outletPos = 3, cellsize = 500)
        ocn_l=landscape_OCN(ocn,slope=0.01)
        ocn_a=aggregate_OCN(ocn_l, thrA = 1.25e6)
        ocn_g=OCN_to_igraph(ocn_a,level = "FD")
        #plot(ocn_g,vertex.size=5)
        m=as.matrix(as_adjacency_matrix(ocn_g))
        g=graph_from_adjacency_matrix(m+t(m),mode="undirected")
        m=adjacency_matrix(g)
        #m=lattice_adjacency_matrix(10,2,1)
        ad=matrix_to_column(m)
        write.table(ad,paste(i,"/matrices/adjacency_matrix_",r, ".txt",sep=""),row.names = F,col.names = F)
  }
