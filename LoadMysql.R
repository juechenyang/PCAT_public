library(RMySQL)
connectDB <- function(){
    options(mysql = list(
      "host" = "private",
      "port" = 0000,
      "user" = "private",
      "password" = "private"
    ))
    databaseName <- "websites"
    db <- dbConnect(MySQL(), dbname = databaseName, host = options()$mysql$host,
      port = options()$mysql$port, user = options()$mysql$user,
      password = options()$mysql$password)
}
# Connect to the database
QueryAll <- function(table) {
  db <- connectDB()
  # Construct the fetching query
  query <- sprintf("SELECT * FROM %s", table)
  # Submit the fetch query and disconnect
  data <- dbGetQuery(db, query)
  dbDisconnect(db)
  return(data)
}
savedata <- function(df, table_name){
  db <- connectDB()
  dbWriteTable(conn = db, name = table_name, value = df, row.names=F, append = TRUE)
}

QueryOnRow <- function(tb_name, row_f){
    db <- connectDB()
    row_f <- paste(row_f,collapse="','")
    query <- paste("select * from ", tb_name, " where gene_name IN ('",row_f,"')", sep="")
    # Submit the fetch query and disconnect
    data <- dbGetQuery(db, query)
    dbDisconnect(db)
    return(data)
}
QueryOnColumn <- function(tb_name, col_f){
    db <- connectDB()
    col_f <- paste(col_f,collapse=",")
    query <- paste("select", col_f, "from", tb_name, sep=" ")
    # Submit the fetch query and disconnect
    data <- dbGetQuery(db, query)
    dbDisconnect(db)
    return(data)
}
DropTable <- function(tb_name){
    db <- connectDB()
    query <- paste("drop table ", tb_name, ";", sep="")
}

