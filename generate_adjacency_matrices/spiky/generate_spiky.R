#generate spiky networks
source("../functions.R")
set.seed(100)
av_degree=c(2) #av. degree of the networks
number_of_nodes=100
number_of_spikes=68 #number of degree 1 patches degree 1 patches
av_degree_ran=(number_of_nodes*av_degree-number_of_spikes)/(number_of_nodes-number_of_spikes) # calculate av. degree of Erdos Renyi network
REPLICATES=1000
for(i in 1:length(av_degree))
{

  system(paste("mkdir ",av_degree[i],sep=""))
  system(paste("mkdir ",av_degree[i],"/matrices",sep=""))
  for(r in 0:(REPLICATES-1))
  repeat
  {
    #generate small Erdos Renyi network
    g<-erdos.renyi.game(number_of_nodes-number_of_spikes,av_degree_ran[i]/(number_of_nodes-number_of_spikes),type="gnp")
    if(is_connected(g)){
      m=adjacency_matrix(g)
      m_spike=matrix(rep(0,(number_of_nodes-number_of_spikes)*number_of_spikes),number_of_spikes)
      spike_pos=sample(1:(number_of_nodes-number_of_spikes),number_of_spikes,replace=T) #pick patches randomly with relacement and add degree 1 patches to them as neighbours
      for(k in 1:length(spike_pos))
      {
        m_spike[k,spike_pos[k]]=1
      }
      m=rbind(m,m_spike)
      m=cbind(m,rbind(t(m_spike),matrix(rep(0,number_of_spikes*number_of_spikes),number_of_spikes)))
      
      mean(apply(m,1,sum))
    
      #m=lattice_adjacency_matrix(10,2,1)
      ad=matrix_to_column(m)
      write.table(ad,paste(av_degree[i],"/matrices/adjacency_matrix_",r, ".txt",sep=""),row.names = F,col.names = F)
      break
    }
  }
}
