#generate maximally modular networks
rm(list=ls()) #clear workspace
source("../functions.R") #load functions required 
#install.packages("igraph") in case it is not installed
library(igraph) # load igraph library
number_of_nodes=96 #number of patches
number_in_module=6 #number of patches in module
#matrix of zeros
m<-matrix(0,nrow=number_of_nodes,ncol=number_of_nodes)
#go through every number_of module-th patch
for(x in seq(1,number_of_nodes-1,number_in_module))
{
  #connect this patch to the number_of module subsequent patches
  m[(x):(x+number_in_module-1),x:(x+number_in_module-1)]=1
  if(x+number_in_module+1<number_of_nodes) 
  {
    #connect one patch to the patch in the next module 
    m[x,(x+number_in_module+1)]=1 
    m[(x+number_in_module+1),x]=1
  }
}

#make a circle
m[2,number_of_nodes]=1
m[number_of_nodes,2]=1

#convert graph to adjacency matrix the function required in in functions.R
ad=matrix_to_column(m)

#write adjacency matrix
write.table(ad,paste(1,"/matrices/adjacency_matrix_",0, ".txt",sep=""),row.names = F,col.names = F)