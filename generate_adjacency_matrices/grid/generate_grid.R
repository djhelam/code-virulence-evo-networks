#generate hexagonal grids
#functions to generate various grids are found in functions.R
source("../functions.R")

neighbours=c(6) #number neighbors of each node
number_of_nodes=100 #number of nodes
n_x=10
n_y=10
for(i in 1:length(neighbours))
{
  m=generate_grids(number_of_nodes,n_x,n_y,neighbours[i] )
  ad=matrix_to_column(m)
  system(paste("mkdir ",neighbours[i],sep=""))
  system(paste("mkdir ",neighbours[i],"/matrices",sep=""))
  write.table(ad,paste(neighbours[i],"/matrices/adjacency_matrix_0.txt",sep=""),row.names = F,col.names = F)
  
}