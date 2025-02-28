  library("fields")
  library("RColorBrewer")
  library("animation")
  library("graphics")
  library("psych")
  #library("modeest")

  # =============================== the rank function ========================================= #
  rank_fun = function(Dmin,Dmax,eta,rank) {
    rank_fun = (Dmax*Dmin*exp(eta*rank)) /  (Dmax+Dmin*(exp(eta*rank)-1))
  }
  r=seq(0,1,by=0.01)
  x.lab=bquote(paste(r[i],"(t)",sep=""))
  y.lab=bquote(paste("D(",r[i],"(t))",sep=""))
  par(mar=c(4,6,1,3))
  plot(r,rank_fun(0.01,19,8,r),type="l",xlab=x.lab,ylab=y.lab,lwd=3,cex.lab=2,cex.axis=2)
  lines(r,rank_fun(0.01,19,5,r),type="l",lwd=3,lty=2)
  lines(r,rank_fun(0.01,2,50,r),type="l",lwd=3,lty=3)
  lines(r,rank_fun(0.01,2,5,r),type="l",lwd=3,lty="431313")
  text(0.12,1.9,bquote(paste(eta,"=",50,sep="")),cex=1.5)
  text(0.26,1.5,bquote(paste(eta,"=",30,sep="")),cex=1.5)
  text(0.58,1,bquote(paste(eta,"=",10,sep="")),cex=1.5)
  text(0.95,0.58,bquote(paste(eta,"=",5,sep="")),cex=1.5)
  # =========================================================================================== #

  # =============================== init ========================================= #
  setwd("~/Documents/scripts/rank")
  #setwd("~/run")

  N=100 #number of agents
  T=500 #number of time steps
  M=1 #number of realizations
  D=0
  eta=0
  #Truth = 1.7
  Truth = 6
  estimates_matrix=matrix(0,N,(T+1))
  #data_unzipped=read.delim("Q-1_Dmax-1_eta-35.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("Q-1_Dmax-1_eta-1.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("Q-1_Dmax-15.5_eta-25.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("P-1_Dmax-7.5_eta-25.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-2_Dmax-10.0_eta-25.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-10.0_eta-25.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-2_Dmax-10.0_eta-10.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-2_Dmax-2.5_eta-5.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-10.0_eta-10.0.rankout.gz",header=FALSE)
  data_unzipped=read.delim("O-1_Dmax-2.5_eta-5.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-20.0_eta-19.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-2.5_eta-19.0.rankout.gz",header=FALSE)
  #data_unzipped=read.delim("O-1_Dmax-14.0_eta-19.0.rankout.gz",header=FALSE)
  # =========================================================================================== #
  
  for (i in seq(1,((T+1)*(N+1)),by=(N+1))) {
    #every (i-1)*(N+1) line is the current time step + additional info on that line
    #everything between (i-1)*(N+1) and i*(N+1) is the estimates at that time step
    estimates_matrix[,(data_unzipped[i,2]+1)] = data_unzipped[(i+1):((i+N)),1]
  }
  ani.options("nmax"=(T+1))

  #plot the histogram and collective error over time
  collective_error=NULL
  Median = NULL
  Mode = NULL
  hist_collective_error.fun=function() {
    for (i in 1:ani.options("nmax")) {
      par(mfrow=c(2,2))  
      hist(estimates_matrix[,i],breaks=70,freq=TRUE,main=paste("Time = ",i,sep=""),xlab="Estimates",ylab="Frequency")
      abline(v=Truth)
      collective_error = c(collective_error,(mean(estimates_matrix[,i])-Truth)^2)      
      plot(1:length(collective_error),collective_error,type="l",col="blue")
      Median = c(Median,(median(estimates_matrix[,i])-Truth)^2)
      plot(1:length(Median),Median,type="l",col="green",ylab="Median",xlab="Time")
      ani.pause(0.001)
    }
  }
 saveMovie(hist_collective_error.fun(),movie.name="test.gif",img.name="Rplot",convert="convert",clean=TRUE)

  #plot the absolute distance between the best agent *at the start* and the truth, and its rank over time
  x11()
  par(mar=c(5,4,4,6)) 
  plot(seq(1:(T+1)),abs(data_unzipped[seq(1,((T+1)*(N+1)),by=(N+1)),5]-Truth),type="l",xlab="Time",ylab="Abs. Distance from Truth",main="Best agent at the start")
  #add the rank of the best agent *at the start* over time. naturally at time 0 rank is 0
  par(new=TRUE)  plot(seq(1:(T+1)),data_unzipped[seq(1,((T+1)*(N+1)),by=(N+1)),6],type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="")
  axis(4)
  mtext("Rank",side=4,line=3)
  
  # plot the absolute distance between the best agent *currently* and the truth over time.
  # *semi-log scale*
  
  for (i in 1:M) {    
    setwd(paste("~/Documents/scripts/rank/Q-",i,sep=""))
    rank_dir = dir(pattern="*rankout*",recursive=TRUE)
      for (j in 1:length(rank_dir)) {
        tokenizer=unlist(strsplit(rank_dir[j],"[\\_,/,-]"))
        Dmax=tokenizer[2];eta=tokenizer[4];m=tokenizer[6];
        data_unzipped=read.delim(gzfile(rank_dir[j]),header=FALSE)
        plot.title=paste("Q-",m,"Dmax-",Dmax,"eta-",eta,sep="")
        pdf(file=paste("../figs/",plot.title,".pdf",sep=""),title=plot.title)
        par(mar=c(5,4,4,6))
          plot(seq(1:(T+1)),log(abs(data_unzipped[seq(1,((T+1)*(N+1)),by=(N+1)),7]-Truth)),type="l",xlab="Time",ylab="Abs. Distance from Truth",main="Best agent currently")
        dev.off()
      }
  }
  
  #add the *mode* and its rank over time. need to bin the data first, else all estimates are different and 'mode' is undefined
  bin_size = 0.1
  mode_distance=NULL
  hist_mode.fun=function() {
     for (i in 1:(T+1)) {
        par(mfrow=c(2,1))
        hist_breaks=seq(min(estimates_matrix[,i]),max(estimates_matrix[,i])+bin_size,by=bin_size)
        h=hist(estimates_matrix[,i],breaks=hist_breaks,freq=TRUE,main=paste("Time = ",i,sep=""),xlab="Estimates",ylab="Frequency")
        h.counts = h$counts
        h.mids = h$mids
        Mode = h.mids[which(h.counts==max(h.counts))]
        abline(v=Mode,col="blue")
        abline(v=Truth)
        mode_distance=c(mode_distance,abs(Mode-Truth))   
        plot(1:length(mode_distance),mode_distance,type="l",col="blue")
        ani.pause(0.001)
    }
  }  
 saveMovie(hist_mode.fun(),movie.name="test.gif",img.name="Rplot",convert="convert",clean=TRUE)
  
  # plot the absolute distance between the mode opinion and the truth over time.
  # *semi-log scale*
   for (i in 1:M) {    
    #setwd(paste("~/Documents/scripts/rank/Q-",i,sep=""))
    setwd(paste("~/run/Q-",1,sep=""))
    rank_dir = dir(pattern="*lnTruth*",recursive=TRUE)
      for (j in 1:length(rank_dir)) {
        tokenizer=unlist(strsplit(rank_dir[j],"[\\_,/,-]"))
        Dmax=tokenizer[2];eta=tokenizer[4];m=tokenizer[6];
        data_unzipped=read.delim(gzfile(rank_dir[j]),header=FALSE)
        for (i in seq(1,((T+1)*(N+1)),by=(N+1))) {
         #every (i-1)*(N+1) line is the current time step + additional info on that line
          #everything between (i-1)*(N+1) and i*(N+1) is the estimates at that time step
          estimates_matrix[,(data_unzipped[i,2]+1)] = data_unzipped[(i+1):((i+N)),1]
        }
        bin_size = 0.1
        mode_distance=matrix(0,1,(T+1))
        for (k in 1:(T+1)) {
          hist_breaks=seq(min(estimates_matrix[,k]),max(estimates_matrix[,k])+bin_size,by=bin_size)
          h=hist(estimates_matrix[,k],breaks=hist_breaks,plot=FALSE)
          h.counts = h$counts
          h.mids = h$mids
          Mode = h.mids[which(h.counts==max(h.counts))]
          mode_distance[k] = min(abs(Mode-Truth))
              #plot(1:length(mode_distance),mode_distance,type="l",col="blue")
        }
        pdf(file=paste("../figs/",plot.title,".pdf",sep=""),title=plot.title)
        plot.title=paste("Q-",1,"Dmax-",Dmax,"eta-",eta,"-mode",sep="")
        par(mar=c(5,4,4,6))
          plot(seq(1:(T+1)),mode_distance,type="l",xlab="Time",ylab="Mode Abs. Distance from Truth",main=plot.title)
        dev.off()
      }
  }

  mode_distance.lagged=diff(mode_distance[1,])
  plot(mode_distance.lagged,type="l")
  skew(mode_distance.lagged)
  mode_distance.lagged.lagged=diff(mode_distance.lagged)
  plot(mode_distance.lagged.lagged,type="l")

  