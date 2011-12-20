library("DBI")
library("RSQLite")
library("fields")
library("RColorBrewer")
library("graphics")

setwd("~/Programs/PySPGfork/PySPG/examples/rank")

#Connect to the database
spg_db = "results_rank.sqlite.copy"
spg_db = "results_rank.sqlite"
SQLite()
drv = dbDriver("SQLite")
con = dbConnect(drv,dbname=spg_db)
dbListTables(con)



#Set of parameters, must match parameters.dat
Truth = seq(1,10,by=2); lTruth = length(Truth)
Realizations = seq(1,100,1); lRealizations = length(Realizations)
Dmax = seq(1,20,by=0.5); lDmax = length(Dmax)
Eta = seq(1,25,by=0.5); lEta = length(Eta)
#auxillary
result_matrix = matrix(NA,length(Eta),length(Dmax))

basefile = "mode-ss-truth-"

for (truthCounter in c(5,7)) {
  truthFlag=FALSE
  pdf(file=paste(basefile,truthCounter,".pdf",sep=""),title="Mode-Ftp")
  for (DmaxCounter in Dmax) {
    #rather stupid way to deal with the spg rounding of the first values of counters
    currentDmax = paste(DmaxCounter,sep="")
    if (as.integer(DmaxCounter)==DmaxCounter && !identical(DmaxCounter,Dmax[1])) 
      currentDmax = paste(DmaxCounter,".0",sep="")
      
    for(EtaCounter in Eta) {
      print(paste(truthCounter,",",DmaxCounter,",",EtaCounter,sep=""))
      currentEta = paste(EtaCounter,sep="")
      if (as.integer(EtaCounter)==EtaCounter && !identical(EtaCounter,Eta[1])) 
        currentEta = paste(EtaCounter,".0",sep="")
    
      spg_db_query=paste("SELECT lnTruth,B,Dmax,eta,mean_fpt,mode_fpt,median_fpt,mean_ss,mode_ss,median_ss FROM values_set,results WHERE values_set.Dmax=",currentDmax," AND values_set.eta=",currentEta," AND values_set.lnTruth=",truthCounter," AND values_set.id=results.values_set_id",sep="")
      spg_db_query_result=NULL
      query=dbSendQuery(con,spg_db_query)
      if (truthCounter != 7) {spg_db_query_result=fetch(query,n=(length(Realizations)-1))}
      if (truthCounter == 7) {spg_db_query_result=fetch(query,n=7)}
      spg_db_result_nrows = nrow(spg_db_query_result)
      if (spg_db_result_nrows > 0 ) {
        truthFlag=TRUE
        data=as.numeric(spg_db_query_result[,"mode_ss"])
        if (length(data[data != -1]) > 0) {
          result_matrix[which(Eta == EtaCounter),which(Dmax == DmaxCounter)] = mean(data[data!=-1])     
        }
      }
      dbClearResult(query)
    }
  }
  if (truthFlag) {
    plot.title=paste("Mode First Passage Time(15%), Truth=",truthCounter,sep="")
    x.lab=bquote(paste(eta,sep=""))
    y.lab=bquote(paste(D[max],sep=""))
    image.plot(Eta,Dmax,result_matrix,xlab=x.lab,ylab=y.lab,cex.lab=2,cex.axis=2)
  }
  dev.off()
}


#dbClearResult(spg_db_query_result)
#spg_db_exceptions = dbGetException(con)
dbDisconnect(con)