cadas30 = function(FC.vector){
  folders = list.dirs('//ussf-prd-isi01/Brigid/cadas3.0', recursive=F)
  select.files = character()
  for (FC in FC.vector){
    current.file = list.files(paste(folders[grep(FC, folders)], '/CADAS_OUTPUT', sep=''), pattern='AllSummary.csv', full.names=TRUE)
    select.files = c(select.files,current.file)
  }
  select.files = lapply(select.files, read.csv, stringsAsFactors=F)
  return(do.call(rbind,select.files))
}