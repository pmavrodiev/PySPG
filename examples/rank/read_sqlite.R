library("DBI")
library("RSQLite")

setwd("/mnt/ethz/home/run")

#Connect to the database
spg_db = "results_rank.sqlite.copy"
SQLite()
drv = dbDriver("SQLite")
con = dbConnect(drv,dbname=spg_db)

#Set of parameters, must match parameters.dat
Truth = seq(1,10,by=2); lTruth = length(Truth)
Realizations = seq(1,100,1); lRealizations = length(Realizations)
Dmax = seq(1,20,by=0.5); lDmax = length(Dmax)
Eta = seq(1,25,by=0.5); lEta = length(Eta)
Results = seq(1,6,by=1); lResults = length(Results)
#Results=c(mean_fpt,mode_fpt,median_fpt,mean_ss,mode_ss,median_ss)

#mean_fpt=1;mode_fpt=1;median_fpt=1;mean_ss=1;mode_ss=1;median_ss=1
#lmean_fpt=length(mean_fpt);lmode_fpt=length(mode_fpt);lmedian_fpt=length(median_fpt);lmean_ss=length(mean_ss);lmode_ss=length(mode_ss);lmedian_ss=length(median_ss)
#Read the spg sqlite db into this multidimensional monster
dim_monster = array(0,dim=c(lTruth,lRealizations,lDmax,lEta,lResults),dimnames=list("Truth"=Truth,"Realizations"=Realizations,"Dmax"=Dmax,"Eta"=Eta,"Results"=Results))


spg_db_query = "SELECT lnTruth,B,Dmax,eta,mean_fpt,mode_fpt,median_fpt,mean_ss,mode_ss,median_ss FROM values_set,results WHERE results.values_set_id = values_set.id"
spg_db_query_result=dbGetQuery(con,spg_db_query)

for (i in 1:nrow(spg_db_query_result)) {
  idx_Truth=which(Truth==spg_db_query_result[i,"lnTruth"])
  idx_Realization=which(Realizations==spg_db_query_result[i,"B"])
  idx_Dmax=which(Dmax==spg_db_query_result[i,"Dmax"])
  idx_Eta=which(Eta==spg_db_query_result[i,"eta"])
  #mean_fpt
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[1]]=spg_db_query_result[i,"mean_fpt"]
  #mode_fpt
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[2]]=spg_db_query_result[i,"mode_fpt"]
  #median_fpt
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[3]]=spg_db_query_result[i,"median_fpt"]
  #mean_ss
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[4]]=spg_db_query_result[i,"mean_ss"]
  #mode_ss
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[5]]=spg_db_query_result[i,"mode_ss"]
  #median_ss
  dim_monster[idx_Truth,idx_Realization,idx_Dmax,idx_Eta,Results[6]]=spg_db_query_result[i,"median_ss"]  
}




dbClearResult(spg_db_query_result)
spg_db_exceptions = dbGetException(con)
dbDisconnect(con)