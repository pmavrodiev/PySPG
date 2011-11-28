library("spam",lib.loc="/home/pmavrodiev/Programs/R/packages")
library("fields",lib.loc="/home/pmavrodiev/Programs/R/packages")
library("RColorBrewer")

  setwd("~/Documents/scripts/rank")
  
  N=100 #number of agents
  T=300 #number of time steps
  M=100 #number of realizations
  W=0
  eta=0
  W_eta_pairs_num=120
  m=0 #current realization
  if (!exists("rank_matrix")) #init only once
    rank_matrix=array(0,dim=c(T,3,W_eta_pairs_num)) #W_eta_pairs {W,eta} pairs
  if (!exists("W_eta_pairs")) 
    W_eta_pairs=matrix(0,W_eta_pairs_num,2)

  #miscellaneous
  DEBUG=1
  IGNORE_ERRORS = 1
  warnings=0
  errors=0

  for (i in 1:M) {

    setwd(paste("~/Documents/scripts/rank/R-",i,sep=""))
    rank_dir = dir(pattern="*rankout*",recursive=TRUE)
    if (length(rank_dir) != W_eta_pairs_num){
      if (DEBUG) {
          print(paste("Error: Number of {W,eta} pairs ",length(rank_dir)," in ",getwd()," does not equal ", W_eta_pairs_num,". Realization skipped",sep=""))
          errors = errors + 1
      }
      next   
    }
   
    for (j in 1:length(rank_dir)) {
      if (exists("DONE_RANK")) break
      #Extract parameters from the file name by tokenizing it using R regular expressions
      #The file name has the following format:
      #W-1.1/eta-0.5/R-1_W-1.1_eta-0.5.rankout.gz
      tokenizer=unlist(strsplit(rank_dir[j],"[\\_,/,-]"))
      W=tokenizer[2];eta=tokenizer[4];m=tokenizer[6];
      #====================Sanity checks====================#
      if (m != i) {
        print("Failed sanity checks")
        errors = errors + 1
        next
      }
      #=====================================================#
      W_eta_pairs[j,1]=W
      W_eta_pairs[j,2]=eta
      #============Process the simulation results===========#
      #a matrix with N rows and 3 columns
      data_unzipped=read.delim(gzfile(rank_dir[j]),header=FALSE)
    

      rank_matrix[,1,j] = rank_matrix[,1,j]+data_unzipped[,1] #later average over all
      rank_matrix[,2,j] = rank_matrix[,2,j]+data_unzipped[,2] #later average over all
      rank_matrix[,3,j] = rank_matrix[,3,j]+data_unzipped[,3] #later average over allrealizations
      #====================Sanity checks====================#
      if (nrow(data_unzipped) != T) {        
          if (!IGNORE_ERRORS) {
            if (DEBUG) {
              print(paste("Error: Number of results ",nrow(data_unzipped)," in ",rank_dir[j]," does not equal T. File skipped.",sep=""))
              errors = errors + 1
            }
          next   
          }
          else {
            if (DEBUG) {
              print(paste("Warning: Number of results ",nrow(data_unzipped)," in ",rank_dir[j]," does not equal T. Discarding unused entries.",sep=""))
              warnings = warnings + 1
            }
          }         
      }
      #=====================================================#
    } #end for (j in 1:length(rank_dir))
   
  } #end for in 1:M

  DONE_RANK = 1

  #average over all realization
  rank_matrix = rank_matrix / M
  setwd("~/Documents/scripts/rank/figs")
  for (k in 1: W_eta_pairs_num) {
    #Plot the collective error and group diversity  
    plot.title=paste("W=",W_eta_pairs[k,1],",eta=",W_eta_pairs[k,2],sep="")  
    pdf(file=paste(plot.title,".pdf",sep=""),title=plot.title)
#    split.screen(c(1,2))
 #   screen(1)

    main.title=plot.title
    plot(rank_matrix[,1,k],rank_matrix[,2,k],type="l",lwd=2,cex.lab=2,cex.axis=2,xlab="t",ylab="Collective Error")
  #  screen(2)
plot(rank_matrix[,1,k],rank_matrix[,3,k],type="l",lwd=2,cex.lab=2,cex.axis=2,xlab="t",ylab="Group Diversity")

   # close.screen(c(1,2))
    dev.off()
    #===================================================#
  }

  #========================================================#
  #=====ALTERNATIVE PLOTTING===============================#
  #PLOT A 2D HEAT MAP WITH time AND eta as x/y axis AND FIXED W#

  t_idx = seq(1,T,by=1)
  eta_idx = seq(0.5,5.0,by=0.5)
  W_eta_pair_matrix = matrix(0,length(t_idx),length(eta_idx))
  W_idx_vector = seq(1.1,2.236,by=0.1)
  setwd("~/Documents/scripts/rank/figs")
  pdf(file=paste("W-const-image-plot.pdf",sep=""),title="W is const")
for (W_idx in seq(1,W_eta_pairs_num,by=10)) {  
    counter = 1
    for (j in W_idx:(W_idx+9))  {
      W_eta_pair_matrix[,counter] = rank_matrix[,2,j] 
      counter = counter + 1
    }
    plot.title = paste("W=",W_idx_vector[(W_idx%/%10 +1)])
image.plot(t_idx,eta_idx,W_eta_pair_matrix,main=plot.title,xlab="Time",ylab="eta")
  }
  dev.off()
    
    
    

  #========================================================#
 