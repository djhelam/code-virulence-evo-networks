rm(list = ls())
#function to generate grids of various neighbour sizes---2,4,6,8
library(igraph) #load packages

generate_grids=function(n,n_x,n_y,neighbours)
{
  m=matrix(rep(0,n*n),n) #stores the adjacency matrix to the graph
  #make a circular graph
  if(neighbours==2)
  {
    x=1:n #patch identities
    for(i in x)#go through all patches
    {
      if(i<n) #two patches are connected if they are neighbour on a line
        m[i,(i+1)]=1 
      if(i>1)
        m[i,(i-1)]=1
      if(i==n)  #periodic boundary conditions
        m[i,1]=1
      if(i==1)
        m[i,n]=1
    }
  }
  #make a grid with 4 nearest neighbours
  if(neighbours==4)
  {
    x=1:n_x #in x-direction
    y=1:n_y  #in y-direction
    all_patches=expand.grid(x=x,y=y) # all x,y coordinates, row number is the node number
    all_patches=cbind(n=1:n,all_patches)
    for(i in 1:n) #go through all rows of all_patches, go through all nodes
    {
      neigh=array() #stores neighbours
      
      #stores the node to the right
      if(all_patches$x[i]<n_x) {neigh[1]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]]
      } else neigh[1]=all_patches$n[all_patches$x==1 & all_patches$y==all_patches$y[i]]
      
      #stores the node to the left
      if(all_patches$x[i]>1){neigh[2]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]]
      } else neigh[2]=all_patches$n[all_patches$x==n_x & all_patches$y==all_patches$y[i]]
      
      #stores the node to the top
      if(all_patches$y[i]<n_y){neigh[3]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]+1]
      } else neigh[3]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==1]
      
      #stores the node to the bottom
      if(all_patches$y[i]>1) {neigh[4]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]-1]
      } else neigh[4]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==n_y]
      
      m[all_patches$n[i],neigh]=1 #if two nodes are connected adjacency matrix element is 1
      
    }
  }
  if(neighbours==8)
  {
    #make a grid with 8 nearest neighbours
    x=1:n_x #in x-direction
    y=1:n_y  #in y-direction
    all_patches=expand.grid(x=x,y=y) # all x,y coordinates, row number is the node number
    all_patches=cbind(n=1:n,all_patches) #node identity
    m=matrix(rep(0,n*n),n) #stores the adjacency matrix to the graph
    for(i in 1:n) #go through all rows of all_patches, go through all nodes
    {
      neigh=array() #stores neighbours
      
      #stores the node to the right
      if(all_patches$x[i]<n_x) {neigh[1]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]]
      } else neigh[1]=all_patches$n[all_patches$x==1 & all_patches$y==all_patches$y[i]]
      
      #stores the node to the left
      if(all_patches$x[i]>1){neigh[2]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]]
      } else neigh[2]=all_patches$n[all_patches$x==n_x & all_patches$y==all_patches$y[i]]
      
      #stores the node to the top
      if(all_patches$y[i]<n_y){neigh[3]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]+1]
      } else neigh[3]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==1]
      
      #stores the node to the bottom
      if(all_patches$y[i]>1) {neigh[4]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]-1]
      } else neigh[4]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==n_y]
      
      #stores the node to the bottom right
      if(all_patches$x[i]<n_x & all_patches$y[i]<n_y) {neigh[5]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]+1] 
      } else {
        if(all_patches$x[i]==n_x & all_patches$y[i]<n_y) {neigh[5]=all_patches$n[all_patches$x==1 & all_patches$y==all_patches$y[i]+1] }
        if(all_patches$x[i]<n_x & all_patches$y[i]==n_y) {neigh[5]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==1] }
        if(all_patches$x[i]==n_x & all_patches$y[i]==n_y) {neigh[5]=all_patches$n[all_patches$x==1 & all_patches$y==1] }
      }
      
      #stores the node to the bottom left
      if(all_patches$x[i]>1 & all_patches$y[i]<n_y) {neigh[6]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]+1]
      } else {
        if(all_patches$x[i]==1 & all_patches$y[i]<n_y) {neigh[6]=all_patches$n[all_patches$x==n_x & all_patches$y==all_patches$y[i]+1] }
        if(all_patches$x[i]>1 & all_patches$y[i]==n_y) {neigh[6]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==1] }
        if(all_patches$x[i]==1 & all_patches$y[i]==n_y) {neigh[6]=all_patches$n[all_patches$x==n_x & all_patches$y==1] }
      }
      
      #stores the node to the top right
      if(all_patches$x[i]<n_x & all_patches$y[i]>1) {neigh[7]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]-1] 
      } else {
        if(all_patches$x[i]<n_x & all_patches$y[i]==1) {neigh[7]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==n_y] }
        if(all_patches$x[i]==n_x & all_patches$y[i]>1) {neigh[7]=all_patches$n[all_patches$x==1 &  all_patches$y==all_patches$y[i]-1]}
        if(all_patches$x[i]==n_x & all_patches$y[i]==1) {neigh[7]=all_patches$n[all_patches$x==1 & all_patches$y==n_y] }
      }
      
      #stores the node to the top left
      if(all_patches$x[i]>1 & all_patches$y[i]>1) {neigh[8]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]-1]
      } else {
        if(all_patches$x[i]>1 & all_patches$y[i]==1) {neigh[8]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==n_y] }
        if(all_patches$x[i]==1 & all_patches$y[i]>1) {neigh[8]=all_patches$n[all_patches$x==n_x &  all_patches$y==all_patches$y[i]-1] }
        if(all_patches$x[i]==1 & all_patches$y[i]==1) {neigh[8]=all_patches$n[all_patches$x==n_x & all_patches$y==n_y] }
      }
      
      m[all_patches$n[i],neigh]=1 #if two nodes are connected adjacency matrix element is 1
      
    }
  }
  if(neighbours==6)
  {
    #make a hexagonal grid with 6 nearest neighbours
    n=100 # number of nodes in the network
    n_x=10 #number of nodes in the x-direction
    n_y=10 #number of nodes in the y-direction
    x_odd=seq(1,n_x,2) #in x-direction
    y_odd=seq(1,2*n_y,2)  #in y-direction
    x_even=seq(2,n_x,2)
    y_even=seq(2,2*n_y,2)
    all_patches=rbind(expand.grid(x=x_odd,y=y_odd),expand.grid(x=x_even,y=y_even)) # all x,y coordinates, row number is the node number
    all_patches=cbind(n=c(seq(1,n,2),seq(2,n,2)),all_patches)
    plot(all_patches$x,all_patches$y)
    m=matrix(rep(0,n*n),n) #stores the adjacency matrix to the graph
    
    for(i in all_patches$n) #go through all rows of all_patches, go through all nodes
    {
      neigh=array() #stores neighbours
      #stores the node to the bottom right
      if(all_patches$x[i]<n_x & all_patches$y[i]<2*n_y) {neigh[1]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]+1]
      } else {
        if(all_patches$x[i]==n_x & all_patches$y[i]<2*n_y) {neigh[1]=all_patches$n[all_patches$x==1 & all_patches$y==all_patches$y[i]+1] }
        if(all_patches$x[i]<n_x & all_patches$y[i]==2*n_y) {neigh[1]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==1] }
        if(all_patches$x[i]==n_x & all_patches$y[i]==2*n_y) {neigh[1]=all_patches$n[all_patches$x==1 & all_patches$y==1] }
      }
      
      #stores the node to the bottom left
      if(all_patches$x[i]>1 & all_patches$y[i]<2*n_y) {neigh[2]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]+1]
      } else {
        if(all_patches$x[i]==1 & all_patches$y[i]<2*n_y) {neigh[2]=all_patches$n[all_patches$x==n_x & all_patches$y==all_patches$y[i]+1] }
        if(all_patches$x[i]>1 & all_patches$y[i]==2*n_y) {neigh[2]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==1] }
        if(all_patches$x[i]==1 & all_patches$y[i]==2*n_y) {neigh[2]=all_patches$n[all_patches$x==n_x & all_patches$y==1] }
      }
      
      #stores the node to the top right
      if(all_patches$x[i]<n_x & all_patches$y[i]>1) {neigh[3]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==all_patches$y[i]-1]
      } else {
        if(all_patches$x[i]<n_x & all_patches$y[i]==1) {neigh[3]=all_patches$n[all_patches$x==all_patches$x[i]+1 & all_patches$y==2*n_y] }
        if(all_patches$x[i]==n_x & all_patches$y[i]>1) {neigh[3]=all_patches$n[all_patches$x==1 &  all_patches$y==all_patches$y[i]-1] }
        if(all_patches$x[i]==n_x & all_patches$y[i]==1) {neigh[3]=all_patches$n[all_patches$x==1 & all_patches$y==2*n_y] }
      }
      
      #stores the node to the top left
      if(all_patches$x[i]>1 & all_patches$y[i]>1) {neigh[4]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==all_patches$y[i]-1]
      } else {
        if(all_patches$x[i]>1 & all_patches$y[i]==1) {neigh[4]=all_patches$n[all_patches$x==all_patches$x[i]-1 & all_patches$y==2*n_y] }
        if(all_patches$x[i]==1 & all_patches$y[i]>1) {neigh[4]=all_patches$n[all_patches$x==n_x &  all_patches$y==all_patches$y[i]-1] }
        if(all_patches$x[i]==1 & all_patches$y[i]==1) {neigh[4]=all_patches$n[all_patches$x==n_x & all_patches$y==2*n_y] }
      }
      
      #stores the node to the bottom
      if( all_patches$y[i]<(2*n_y-1)) {neigh[5]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]+2]
      }  else {
        if(all_patches$y[i]==(2*n_y)){neigh[5]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==(2)]}
        if(all_patches$y[i]==(2*n_y-1)){neigh[5]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==(1)]}
      }
      
      #stores the node to the top
      if( all_patches$y[i]>2) {neigh[6]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==all_patches$y[i]-2]
      } else {
        if(all_patches$y[i]==1){neigh[6]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==(2*n_y-1)]}
        if(all_patches$y[i]==2){neigh[6]=all_patches$n[all_patches$x==all_patches$x[i] & all_patches$y==(2*n_y)]}
      }
      m[all_patches$n[i],neigh]=1 #if two nodes are connected adjacency matrix element is 1
    }
  }
  
  return(m)
  
  
}


generate_rgg=function(n, radius)
{
  g=sample_grg(n,radius,coor=T) #RGG with n nodes and specified radius
  return(g)
}

adjacency_matrix=function(g)   #create a random geometric graph RGG with a specified radius
{
  
  return(as.matrix(as_adjacency_matrix(g))) #return adjacency matrix
}


matrix_to_column=function(mat) #each element of the adjacency matrix in a new line
{
  column=data.frame()
  for(i in 1:nrow(mat))
  {
    temp=data.frame(mat[i,])
    column=rbind(column,temp)
  }
  return(column)
}

