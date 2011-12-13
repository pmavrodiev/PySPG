  library("fields")
  library("RColorBrewer")
  library("animation")
  library("graphics")
  
  # ========================== INITIALIZATION ===========================#
  root_dir_name="~/Documents/scripts/rank/"
  setwd(root_dir_name)
  N=100 #number of agents
  T=500 #number of time steps
  M=9 #number of realizations
  eta=0
  Truth = -0.125
  estimates_matrix=matrix(0,N,(T+1))
  root_dir_name = "~/run/Q-"
  bin_size = 0.1
  Truth_bracket = abs(0.15*Truth)
  Steady_state_bracket=0.1
  Skew_boundary=1 #convergence boundary; if skew larger than the boundary then the series converge
  # ===================================================================#
  #file_name = "Q-1_Dmax-1_eta-5.0.rankout.gz"
  for (i in 1:M) {
    print(paste("Realization ",i,sep=""))
    setwd(paste(root_dir_name,"/Q-",i,sep=""))
    # === Extract the number of Dmax values and the values themselves ===#
    rank_dir = dir(pattern="*Dmax*",recursive=FALSE)
    Dmax_sequence = matrix(0,1,length(rank_dir))
    for (j in 1:length(rank_dir)) {
      tokenizer = unlist(strsplit(rank_dir[j],"[-]"))
      Dmax_sequence[j]=as.numeric(tokenizer[2])
    }
    Dmax_sequence = sort(Dmax_sequence)
    # ====== Create a 'list' for faster retrieval of indeces later ======#
    Dmax_index_list = list() # empty list
    for (k in 1:length(Dmax_sequence)) 
      Dmax_index_list[[as.character(Dmax_sequence[k])]] = k 
    # ===================================================================#
    # === Extract the number of eta values and the values themselves ===# 
    setwd(paste(root_dir_name,"/Q-",i,"/Dmax-1",sep=""))
    rank_dir = dir(pattern="*eta*",recursive=FALSE)
    eta_sequence = matrix(0,1,length(rank_dir))
    for (j in 1:length(rank_dir)) {
      tokenizer = unlist(strsplit(rank_dir[j],"[-]"))
      eta_sequence[j]=as.numeric(tokenizer[2])
    }
    eta_sequence = sort(eta_sequence)
    # ====== Create a 'list' for faster retrieval of indeces later ======#
    eta_index_list = list() # empty list
    for (k in 1:length(eta_sequence)) 
      eta_index_list[[as.character(eta_sequence[k])]] = k
    # ===================================================================#
    Dmax_eta_matrix = matrix(NA,length(eta_sequence),length(Dmax_sequence))
    #number of realizations for each Dmax,eta pair for which data can be computed 
    Dmax_eta_matrix_number = matrix(0,length(eta_sequence),length(Dmax_sequence))
    # ==================== start reading files ==========================#
    MEDIAN_DISTANCE=TRUE
    MODE_DISTANCE=FALSE
    STATIONARY_STATE_CONVERGENCE=FALSE
    setwd(paste(root_dir_name,"/Q-",i,sep=""))
    rank_dir = dir(pattern="*rankout*",recursive=TRUE)
    for (j in 1:length(rank_dir)) {
      print(rank_dir[j],sep="")
      data_unzipped=read.delim(rank_dir[j],header=FALSE)
      if (length(data_unzipped[,1]) != 50601 ) next  #101101
      tokenizer=unlist(strsplit(rank_dir[j],"[\\_,/,-]"))
      Dmax=tokenizer[2];eta=tokenizer[4];
      # ======= Read in the individual estimates for all time steps =====#
      for (ii in seq(1,((T+1)*(N+1)),by=(N+1))) 
        estimates_matrix[,(data_unzipped[ii,2]+1)] = data_unzipped[(ii+1):((ii+N)),1]
      # ====================== Calculate the mode =======================#
      if (MEDIAN_DISTANCE) {
        median_distance=matrix(0,1,(T+1))
        for (k in 1:(T+1)) {
          median_distance[k] = median(estimates_matrix[,k])
          if (median_distance[k] <= Truth_bracket) {
            print(paste("Found: Dmax=",Dmax," eta=",eta,sep=""))
            eta = as.character(as.numeric(eta)) # to lose the most significant digit of an int, i.e. from 3.0 to 3
            Dmax = as.character(as.numeric(Dmax))
            oldValue=Dmax_eta_matrix[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]
            newValue=k
            if (!is.na(oldValue)) newValue = oldValue+k
            #save the time when the mode falls within the defined bracket
            Dmax_eta_matrix[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]] = newValue
            Dmax_eta_matrix_number[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]=Dmax_eta_matrix_number[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]+1
            break          
          }
        }
      }
      if (MODE_DISTANCE) {
        mode_distance=matrix(0,1,(T+1))
        for (k in 1:(T+1)) {
          hist_breaks=seq(min(estimates_matrix[,k]),max(estimates_matrix[,k])+bin_size,by=bin_size)
          h=hist(estimates_matrix[,k],breaks=hist_breaks,plot=FALSE)
          h.counts = h$counts
          h.mids = h$mids
          Mode = h.mids[which(h.counts==max(h.counts))]
          mode_distance[k] = min(abs(Mode-Truth)) #in case of multi-modal scenario, take the closest to the Truth
      
            if (mode_distance[k] <= Truth_bracket) {
              print(paste("Found: Dmax=",Dmax," eta=",eta,sep=""))
              eta = as.character(as.numeric(eta)) # to lose the most significant digit of an int, i.e. from 3.0 to 3
              Dmax = as.character(as.numeric(Dmax))
              oldValue=Dmax_eta_matrix[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]
              newValue=k
              if (!is.na(oldValue)) newValue = oldValue+k
              #save the time when the mode falls within the defined bracket
              Dmax_eta_matrix[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]] = newValue
            Dmax_eta_matrix_number[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]=Dmax_eta_matrix_number[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]]+1
              break
            }
        }
      }
      if (STATIONARY_STATE_CONVERGENCE) {
        mode_distance_lagged = diff(mode_distance[1,])
        if (abs(skew(mode_distance_lagged)) >= Skew_boundary) { #convergence
          for (k in 1:(T+1)) {
            steady_state = mode_distance[T+1]
            if (steady_state == 0) print("Pathological case")
            if (abs(mode_distance[k]) <= Steady_state_bracket*steady_state) {
              print(paste("Found2: Dmax=",Dmax," eta=",eta,sep=""))
              eta = as.character(as.numeric(eta)) # to lose the most significant digit of an int, i.e. from 3.0 to 3
              Dmax = as.character(as.numeric(Dmax))
              Dmax_eta_matrix[unlist(eta_index_list[eta])[[1]],unlist(Dmax_index_list[Dmax])[[1]]] = k
              break
            }              
          } 
        }
      }     
      # ===================================================================#
    }
    # ===================================================================#
  } #end for (i in 1:M)
  if (STATIONARY_STATE_CONVERGENCE) {
    plot.title=paste("Mode, Time to Steady State (10%), Truth=",6,sep="")
    image.plot(eta_sequence,Dmax_sequence,Dmax_eta_matrix,main=plot.title,xlab="Sensitivity (eta)",ylab="Maximum diffusion (Dmax)")
  }
  if (MODE_DISTANCE) {
    #average
    Dmax_eta_matrix = Dmax_eta_matrix / Dmax_eta_matrix_number
    plot.title=paste("Mode First Passage Time(15%), Truth=",Truth,sep="")
    
    image.plot(eta_sequence,Dmax_sequence,Dmax_eta_matrix,main=plot.title,xlab="Sensitivity (eta)",ylab="Maximum diffusion (Dmax)")
  }