#generate random networks
source("../functions.R")
set.seed(100) #set seed for random number generator
av_degree=6 # average degree of Erdos Renyi network
number_of_nodes=100 #number of nodes
REPLICATES=1000 #number of realisatiosn
for(i in 1:length(av_degree))
{

  system(paste("mkdir ",av_degree[i],sep=""))
  system(paste("mkdir ",av_degree[i],"/matrices",sep=""))
  for(r in 0:(REPLICATES-1))
  repeat
  {
    #generate Erdos renyi network using igraph function
    g<-erdos.renyi.game(number_of_nodes,av_degree[i]/number_of_nodes,type="gnp")
   #if the graph is connected write the file containing the adjacency matrix
     if(is_connected(g)){
      m=adjacency_matrix(g)
      mean(apply(m,1,sum))
      #m=lattice_adjacency_matrix(10,2,1)
      ad=matrix_to_column(m)
      write.table(ad,paste(av_degree[i],"/matrices/adjacency_matrix_",r, ".txt",sep=""),row.names = F,col.names = F)
      break
    }
  }
}