library("DBI")
library("RSQLite")

setwd("~/run")
SQLite()
drv = dbDriver("SQLite")
con = dbConnect(drv,dbname="results_rank.sqlite")

#list all tables
tab =dbListTables(con)
data = dbReadTable(conn=con,name=tab[3])


dbListTables(con)
dbDisconnect(con)