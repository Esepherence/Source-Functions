#######################################
#FC = flowcell id
#runfolders = where to search for FC
#skip '21' will not work for all sample sheets!!!

samplesheet = function(FC, runfolders){
  current.file = list.files(paste0(runfolders[grep(FC, runfolders)], '/'), pattern='SampleSheet.csv$', full.names=TRUE)
  
  df.temp = read.table(as.character(current.file),header=F)
  skip.to.data = grep("DATA",df.temp$V1,ignore.case=T)
  
  current.file = read.csv(current.file, stringsAsFactors=F, skip=skip.to.data)
  return(current.file)
}

#######################################
#inputs a vector of FCs
#calls the samplesheet function for sampleIDs
#merges sample IDs found in sample sheets with cadas output and returns output

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

#######################################
#takes XX/XY seq data input (most likely form cadas30 function)
#evaluates failing conditions
#returns a list of failed pools

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
#######################################
#Reads google doc

getBrigid = function(){
  ut = require('googlesheets')
  if(ut==F){
    install.packages('devtools')
    devtools::install_github("jennybc/googlesheets")
  }else{
    library('googlesheets')
  }
  gs_auth()
  br <- gs_key(gs_ls()$sheet_key[1])
  getData = function(x){
    df = gs_read_csv(br,ws=x)
    df 
  }
  dd = lapply(as.list(br$ws$ws_title),FUN=getData)
  out = lapply(dd,as.data.frame,header=T)
  names(out) = br$ws$ws_title
  
  return(out)
}

#######################################
#updates parsed SAV files

update.sav = function(){
  require(plyr)
  source('c:/pe_repo/r_pe/tech_transfer/Brigid/googlesBrigid.R')
  br = getBrigid()
  
  #set home settings
  #dirP = '//srv2isilon/hiseq_output/'
  dirP = '//ussf-prd-isi01/brigid/raw_data/'
  outD = '//ussf-prd-isi01/Brigid/rnd/nextseq/'

  #get run folders
  run_paths = function(dirP){
    inst = dir(dirP)[grepl('nextseq[0-9]',dir(dirP))]
    getFC = function(inst){
      path = paste(dirP,inst,sep='')
      f = dir(path)
      f = f[grepl('_A',f)]
      if(length(f)==0) return()
      out = data.frame('Path' = as.character(f),"Inst"=inst,dir = dirP)
      out
    }
    
    dd = lapply(inst,getFC)
    ##remove empty FC data
    dd[sapply(dd,is.null)] = NULL
    # build data frame
    fcl = do.call(rbind.fill,dd)
    fcl
  }

  fcl = run_paths(dirP)
  #add SD data
  dirsd = '//sd-isilon/Brigid/raw_data/'
  fcl2 = run_paths(dirsd)
  fcl = rbind(fcl,fcl2)


  fc = function(x){ strsplit(x,'_')[[1]][4]}
  fcl$fcid = apply(fcl,1,fc)
  fcl$DATE = substr(fcl$Path,1,6)
  fcl$RDATE = strptime(fcl$DATE,forma='%y%m%d')
  fcl$Flowcell = substr(fcl$fcid,2,10)
  mk = fcl$Flowcell %in% br$Sequencing$FLOW_CELL_ID

  runs = fcl[mk,]
  setwd(outD)
  for(i in 1:dim(runs)[1]){
    fName = paste0(runs$fcid[i],'.xml')
    file_out = paste('\"',outD,fName,'\"',sep='')
    #file_outT = paste('\"',outD,runs$fcid[i],'_table','.xml\"',sep='')
    if(!fName %in% dir(outD)){
      sav = '\"C:\\Illumina\\Illumina Sequencing Analysis Viewer\\Sequencing Analysis Viewer.exe\"'
      st = try(system(paste(sav,paste(runs$dir[i],runs$Inst[i],'/',runs$Path[i],sep=''),file_out,"-s")))
      #st = try(system(paste(sav,paste(dirP,runs$Inst[i],'/',runs$Path[i],sep=''),file_outT,"-t")))
      print(paste(sav,paste(runs$dir,runs$Inst[i],'/',runs$Path[i],sep=''),file_out,"-s"))
      print(i)
    }	
  }
}


#######################################
#Takes FC id(s) and returns intensity df

#example: require("plyr");int = ldply("HC5CNBGXX",getInt1)
getInt1 = function(FC){
  setwd("//ussf-prd-isi01/Brigid/rnd/nextseq/")
  filePath = grep(FC,dir(),value=T)
  int= xmlTreeParse(filePath)
  xmltop = xmlRoot(int)
  data = xmlGetAttr(xmltop,'FirstCycleInt')
  q3 = xmlGetAttr(xmltop,'PrcGtQ30')
  out = data.frame(Flowcell=substr(filePath,2,10),C1=data,Q30=q3)
  out
}

#######################################
#Takes FC id(s) and returns lane summary

#example: require("plyr");lan = ldply("HC5CNBGXX",getLane)
# lan$ClustersPF = ave(as.numeric(as.character(lan$ClustersCountPf.Mean))*as.numeric(as.character(lan$.attrs.TileCount)),paste(lan$Flowcell,lan$Read),FUN=sum)
# lan$Density = ave(as.numeric(as.character(lan$Density.Mean))/1000,paste(lan$Flowcell,lan$Read),FUN=mean)
getLane = function(FC){
  setwd("//ussf-prd-isi01/Brigid/rnd/nextseq/")
  filePath = grep(FC,dir(),value=T)
  print(filePath)
  int = xmlParse(filePath)
  ll = xmlToList(int)
  read1 = data.frame(rbind(unlist(ll[['SummaryByRead']][['ReadSummary']][['SummaryByLane']][[1]]),unlist(ll[['SummaryByRead']][['ReadSummary']][['SummaryByLane']][[2]]),unlist(ll[['SummaryByRead']][['ReadSummary']][['SummaryByLane']][[3]]),unlist(ll[['SummaryByRead']][['ReadSummary']][['SummaryByLane']][[4]])))
  if(is.null(dim(read1)) | dim(read1)[1]==0)return()
  if(read1$Density.Mean[1]==0) return()
  read1$Read = 1
  read2 = data.frame(rbind(unlist(ll[['SummaryByRead']][[2]][['SummaryByLane']][[1]]),unlist(ll[['SummaryByRead']][[2]][['SummaryByLane']][[2]]),unlist(ll[['SummaryByRead']][[2]][['SummaryByLane']][[3]]),unlist(ll[['SummaryByRead']][[2]][['SummaryByLane']][[4]])))
  read2$Read = 2
  reads = rbind(read1,read2)
  if(dim(reads)[1]==0) return(reads)
  reads$Flowcell = substr(filePath,2,10)
  reads
}

