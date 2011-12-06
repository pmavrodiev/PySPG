library("fields")
library("RColorBrewer")
library("animation")
library("graphics")
  setwd("~/Documents/scripts/rank")

  N=100 #number of agents
  T=500 #number of time steps
  M=3 #number of realizations
  D=0
  eta=0
  Truth = 1.7
  #Warning: magic numbers here, sequences below must match parameters_rank.dat
  #####################################
  estimates_matrix=matrix(0,N,(T+1))
  collective_error=NULL
  
  data_unzipped=read.delim("O-1_Dmax-10.0_eta-10.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-2.5_eta-5.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-20.0_eta-19.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-2.5_eta-19.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-14.0_eta-19.0.rankout.gz",header=FALSE)
  
  for (i in seq(1,((T+1)*(N+1)),by=(N+1))) {
    #every (i-1)*(N+1) line is the current time step + additional info on that line
    #everything between (i-1)*(N+1) and i*(N+1) is the estimates at that time step
    estimates_matrix[,(data_unzipped[i,2]+1)] = data_unzipped[(i+1):((i+N)),1]
  }
  ani.options("nmax"=(T+1))
  ani.fun=function() {
    for (i in 1:ani.options("nmax")) {
      par(mfrow=c(2,1))  
      hist(estimates_matrix[,i],breaks=70,freq=TRUE,main=paste("Time = ",i,sep=""),xlab="Estimates",ylab="Frequency")
      abline(v=Truth)
      collective_error = c(collective_error,(mean(estimates_matrix[,i])-Truth)^2)      
      plot(1:length(collective_error),collective_error,type="l",col="blue")
      ani.pause(0.001)
    }
  }

  saveMovie(ani.fun(),movie.name="test.gif",img.name="Rplot",convert="convert",clean=TRUE)

    