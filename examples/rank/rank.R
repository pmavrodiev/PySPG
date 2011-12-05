library("fields")
library("RColorBrewer")

  setwd("~/Documents/scripts/rank")
  
  N=100 #number of agents
  T=3000 #number of time steps
  M=3 #number of realizations
  D=0
  eta=0
  t_idx = seq(1,T,by=1)
  #Warning: magic numbers here, sequences below must match parameters_rank.dat
  eta_idx = seq(1,35,by=0.5)
  Dmax_idx_vector = seq(1,20,by=0.5)
  #####################################
  Dmax_eta_pairs_num=length(eta_idx)*length(Dmax_idx_vector)
  Dmax_eta_pair_matrix = matrix(0,length(t_idx),length(eta_idx))
  

  m=0 #current realization
  if (!exists("rank_matrix")) #init only once
    rank_matrix=array(0,dim=c(T,3,Dmax_eta_pairs_num)) #Dmax_eta_pairs {Dmax,eta} pairs
  if (!exists("Dmax_eta_pairs")) 
    Dmax_eta_pairs=matrix(0,Dmax_eta_pairs_num,2)

  #miscellaneous
  DEBUG=1
  IGNORE_ERRORS = 1
  warnings=0
  errors=0

  for (i in 1:M) {

    setwd(paste("~/Documents/scripts/rank/O-",i,sep=""))
    rank_dir = dir(pattern="*rankout*",recursive=TRUE)
    #SANITY CHECKS
    if (length(rank_dir) != Dmax_eta_pairs_num){
      if (DEBUG) {
          print(paste("Error: Number of {Dmax,eta} pairs ",length(rank_dir)," in ",getwd()," does not equal ", Dmax_eta_pairs_num,". Realization skipped",sep=""))
          errors = errors + 1
      }
      PROMPT = ""
      while (PROMPT != "y" && PROMPT != "n")
        PROMPT = readline("Continue? (y/n) ")
      if (substr(PROMPT,1,1) == "y")
        next   
      else break
    }
    ###########################
    for (j in 1:length(rank_dir)) {
      if (exists("DONE_RANK")) break
      #Extract parameters from the file name by tokenizing it using R regular expressions
      #The file name has the following format:
      #W-1.1/eta-0.5/R-1_W-1.1_eta-0.5.rankout.gz
      tokenizer=unlist(strsplit(rank_dir[j],"[\\_,/,-]"))
      
      Dmax=tokenizer[2];eta=tokenizer[4];m=tokenizer[6];
      #====================Sanity checks====================#
      if (m != i) {
        print("Failed sanity checks. This should never happen")
        errors = errors + 1
        next
      }
      #=====================================================#
      Dmax_eta_pairs[j,1]=Dmax
      Dmax_eta_pairs[j,2]=eta
      #============Process the simulation results===========#
      #a matrix with N rows and 3 columns
      data_unzipped=read.delim(gzfile(rank_dir[j]),header=FALSE)
      #====================Sanity checks====================#
      flag=FALSE
      if (nrow(data_unzipped) != T) {        
          if (nrow(data_unzipped) == (2*T)) {
            if (DEBUG) {
              print(paste("Warning: Number of results ",nrow(data_unzipped)," in ",rank_dir[j]," equals 2*T; two spg workers operated on the same file.",sep=""))
              warnings = warnings + 1
              flag=TRUE
            }
          }
          else {
            if (!IGNORE_ERRORS) {
              if (DEBUG) {
                print(paste("Error: Number of results ",nrow(data_unzipped)," in ",rank_dir[j]," does not equal T. File skipped.",sep=""))
                errors = errors + 1
              }
              PROMPT = ""
              while (PROMPT != "y" && PROMPT != "n")
                PROMPT = readline("Continue? (y/n) ")
              if (substr(PROMPT,1,1) == "y")
                next   
              else break
            }
            else {
              if (DEBUG) {
              print(paste("Warning: Number of results ",nrow(data_unzipped)," in ",rank_dir[j]," does not equal T. Discarding file.",sep=""))
                warnings = warnings + 1
              }
            }         
        }
      }
      #=====================================================#
      if (!flag) {
        rank_matrix[1:T,1,j] = rank_matrix[1:T,1,j]+data_unzipped[1:T,1] #later average over all
        rank_matrix[1:T,2,j] = rank_matrix[1:T,2,j]+data_unzipped[1:T,2] #later average over all
        rank_matrix[1:T,3,j] = rank_matrix[1:T,3,j]+data_unzipped[1:T,3] #later average over allrealizations
      }
      else {
        #average the two realizations
        for (ii in 1:T) {
          data_unzipped[ii,2]= (data_unzipped[ii,2]+data_unzipped[ii+T,2])/2
          data_unzipped[ii,3]= (data_unzipped[ii,3]+data_unzipped[ii+T,3])/2
        }       
        rank_matrix[1:T,1,j] = rank_matrix[1:T,1,j]+data_unzipped[1:T,1] #later average over all
        rank_matrix[1:T,2,j] = rank_matrix[1:T,2,j]+data_unzipped[1:T,2] #later average over all
        rank_matrix[1:T,3,j] = rank_matrix[1:T,3,j]+data_unzipped[1:T,3] #later average over allrealizations
      }
    } #end for (j in 1:length(rank_dir))
   
  } #end for in 1:M

  DONE_RANK = 1

  #average over all realization
  rank_matrix = rank_matrix / M
  setwd("~/Documents/scripts/rank/figs")
  for (k in 1: Dmax_eta_pairs_num) {
    #Plot the collective error and group diversity  
    plot.title=paste("Dmax=",Dmax_eta_pairs[k,1],",eta=",Dmax_eta_pairs[k,2],sep="")  
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


  setwd("~/Documents/scripts/rank/figs")
  pdf(file=paste("Dmax-const-image-plot.pdf",sep=""),title="Dmax is const")
for (Dmax_idx in seq(1,Dmax_eta_pairs_num,by=length(eta_idx))) {  
    counter = 1
    for (j in Dmax_idx:(Dmax_idx+length(eta_idx)-1))  {
      Dmax_eta_pair_matrix[,counter] = rank_matrix[,2,j] 
      counter = counter + 1
    }
    plot.title = paste("Dmax=",Dmax_idx_vector[(Dmax_idx%/%length(eta_idx) +1)])
image.plot(t_idx,eta_idx,Dmax_eta_pair_matrix,main=plot.title,xlab="Time",ylab="eta")
  }
  dev.off()
    
    #========================================================#
 