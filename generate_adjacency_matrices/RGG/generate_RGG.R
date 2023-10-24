#this is the R script to generate the 1000 RGGs used in the paper
source("../functions.R") #functions required
set.seed(100) #set seed for random number generator
rad=c(0.15) #radius parameter of the RGG
number_of_nodes=100 #nuber of nodes in the graph
REPLICATES=1000 #number of graphs to generate (realisations)
for(i in 1:length(rad))
{

  system(paste("mkdir ",i,sep="")) #create a directory corresponding to a radius
  system(paste("mkdir ",i,"/matrices",sep="")) #store all adjacency matrices
  for(r in 0:(REPLICATES-1)) #go through all replicates
  #we only want netowrks that are connected
  #so generate RGGs, discard if they are not connected
  repeat 
  {
    #use igraph function to generate the RGG
    g=generate_rgg(number_of_nodes,rad[i])
    #write to folder if the network is connected
    if(is_connected(g)){
      m=adjacency_matrix(g)
      ad=matrix_to_column(m)
      write.table(ad,paste(i,"/matrices/adjacency_matrix_",r, ".txt",sep=""),row.names = F,col.names = F)
      break
    }
  }
}