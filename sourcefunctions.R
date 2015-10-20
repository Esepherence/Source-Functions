samplesheet = function(FC, runfolders){
  current.file = list.files(paste(runfolders[grep(FC, runfolders)], '/', sep=''), pattern='SampleSheet.csv$', full.names=TRUE)
  current.file = read.csv(current.file, stringsAsFactors=F, skip=21)
  return(current.file)
}

cadas30 = function(FC.vector){
  require(plyr)
  folders = list.dirs('//ussf-prd-isi01/Brigid/cadas3.0', recursive=F)
  rawfolders = list.dirs('//ussf-prd-isi01/Brigid/raw_data', recursive=F)
  rawfolders = rawfolders[grepl('next', rawfolders)]
  runfolders = c()
  for (folder in rawfolders){
    runfolders = append(runfolders, list.dirs(folder, recursive=F))
  }
  select.files = data.frame()
  for (FC in FC.vector){
    current.file = list.files(paste(folders[grep(FC, folders)], '/CADAS_OUTPUT', sep=''), pattern='_AllSummary.csv', full.names=TRUE)
    current.file = read.csv(current.file, stringsAsFactors=F)
    samplesheet = samplesheet(FC, runfolders)  
    select.files = rbind.fill(select.files,merge(current.file, samplesheet, by='SampleID', all.x=T, suffixes=c('', 'SS')))
    print(paste('Flowcell Complete: ', FC))
  }
  return(select.files)
}




poolqual = function(data){
  data$PoolSex = substr(data$SampleID, 1, 2)
  poolmk = (data$PoolSex == 'XX' | data$PoolSex == 'XY')
  data = data[poolmk,]
  nesmk = data$NonExcludedSites < 4000000
  mk13 = data$NCV_13 > 2.5
  mk18 = data$NCV_18 > 2.5
  mk21 = data$NCV_21 > 2.5
  mksex = (data$PoolSex == 'XX' & (data$NCV_Y > 5 | data$NCV_Y < -5)) | ((data$PoolSex == 'XY'& data$NCV_Y < 30))
  ffmk = data$FF < .08
  failnes = data$SampleID[nesmk]
  fail13 = data$SampleID[mk13]
  fail18 = data$SampleID[mk18]
  fail21 = data$SampleID[mk21]
  failsex = data$SampleID[mksex]
  failff = data$SampleID[ffmk]
  failed = data$SampleID[(nesmk | mk13 | mk18 | mk21 | mksex | ffmk)]
  return(list(failnes = failnes, fail13 = fail13, fail18 = fail18, fail21 = fail21, failsex = failsex, failff = failff, FAILED = failed))
}