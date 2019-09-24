if(!"BiocManager" %in% rownames(installed.packages())){
  install.packages("BiocManager")
}
library(BiocManager)

packageF <- function(package){
  stat = require(package, character.only = T)
  
  warnReal <- warnings()
  
  if (stat) {
    print("All is good")
  } else {
    BiocManager::install(package)
    require(package, character.only = T)
  }
}

packageF("data.table")
packageF("scales")
packageF("memoise")
packageF("sva")
packageF("lme4")
packageF("car")
packageF("corpcor")
if(!"igraph" %in% installed.packages()){
  install.packages("http://cran.r-project.org/src/contrib/Archive/igraph/igraph_1.2.1.tar.gz", repos=NULL, type="source")
}

packageF("igraph")
packageF("magrittr")
packageF("sm")
packageF("ggplot2")
packageF("dplyr")
packageF("curl")

packageF("gplots")
packageF("scales")
packageF("devtools")
packageF("magrittr")
select = dplyr::select
filter = dplyr::filter
mutate = dplyr::mutate

name2Char <- function(x){
  as.character(substitute(x))
}


GeneSex <- function(aned, Metadata){
  probeM <- aned$Probe[grep("RPS4Y1|KDM5D", aned$GeneSymbol, ignore.case = T)] %>% as.character
  probeF <- aned$Probe[grep("XIST", aned$GeneSymbol, ignore.case = T)] %>% as.character
  
  if(length(probeM) > 0 | length(probeF) > 0){
    dataSex <- t(aned %>%
                   filter(Probe %in% c(probeM, probeF)) %>%
                   dplyr::select(matches("(_)|(GSM)")))
    colnames(dataSex) <- aned$Probe[aned$Probe %in% c(probeM, probeF)]
    Clusters <- kmeans(dataSex, centers=2)
    
    Centers <- Clusters$center
    
    if(length(probeM)>0 & length(probeF)>0){
      if(mean(Centers[1,probeM]) > mean(Centers[1,probeF]) & mean(Centers[2,probeM]) < mean(Centers[2,probeF])){
        Clusters$cluster[Clusters$cluster==1] <- "M"
        Clusters$cluster[Clusters$cluster==2] <- "F"
      } else if(mean(Centers[1,probeM]) < mean(Centers[1,probeF]) & mean(Centers[2,probeM]) > mean(Centers[2,probeF])){
        Clusters$cluster[Clusters$cluster==1] <- "F"
        Clusters$cluster[Clusters$cluster==2] <- "M"
      } else {
        warning("Gender genes disagree, cannot decide about biological gender")
      }
    } else if(length(probeF) == 0) {
      if (mean(Centers[1,]) > mean(Centers[2,])) {
        Clusters$cluster[Clusters$cluster==1] <- "M"
        Clusters$cluster[Clusters$cluster==2] <- "F"
      } else {
        Clusters$cluster[Clusters$cluster==1] <- "F"
        Clusters$cluster[Clusters$cluster==2] <- "M"
      }  
    } else if (length(probeM) == 0) {
      if (mean(Centers[1,]) > mean(Centers[2,])){
        Clusters$cluster[Clusters$cluster==1] <- "F"
        Clusters$cluster[Clusters$cluster==2] <- "M"
      } else {
        Clusters$cluster[Clusters$cluster==1] <- "M"
        Clusters$cluster[Clusters$cluster==2] <- "F"
      }
    } 
    if("Gender genes disagree, cannot decide about biological gender" %in% names(warnings())){
      Metadata$Biogender <- NA
    } else {
      Metadata$BioGender <- Clusters$cluster[match(Metadata$CommonName, names(Clusters$cluster))]
    }
    
    return(Metadata)
  } else {
    warning("No sex specific genes on the platform")
  }
}

closeDev <- function(){
  Dev2close <- dev.list()[names(dev.list()) != "RStudioGD"]
  sapply(Dev2close, function(x) dev.off(x))
}

plm_data_function <- function(rawData){
  PLM_data <- fitProbeLevelModel(temp) #takes a while... can alsa be used on the rma-preproccessed data, but requires to change input
  return(PLM_data)
}


F2N <- function(x){
  return(as.numeric(as.character(x)))
}

GetSigChar <- function(pVal){
  if (pVal > 0.05) {
    sigChar = ""
  } else if (pVal <= 0.05 & pVal > 0.01){
    sigChar = "*"
  } else if (pVal <= 0.01 & pVal > 0.001){
    sigChar = "**"
  } else if (pVal <= 0.001) {
    sigChar ="***"
    return(sigChar)
  }
}

GetSigCol <- function(pVal){
  if (pVal < 0.05) {
    col = "red"
  } else {
    col = "black"
  }
}

presentSig <- function(pVal, signifChar = TRUE, signifCol = TRUE){
  #Choosing what data to present
  
  if(signifChar & pVal < 0.05) {
    sigChar <- GetSigChar(pVal)
  } else {
    sigChar = ""
  }
  
  if(signifCol){
    sigCol <- GetSigCol(pVal)
  } else {
    sigCol = "black"
  }
  
  if(pVal > 0.05) {
    pVal = "n.s."
  }
  return(list("sigChar" = sigChar,
              "sigCol" = sigCol,
              "pVal" = pVal))
}

List2df <- function(List, sep="\\."){
  df <- do.call(cbind, List) %>% as.data.frame
  if(sum(grepl(sep, names(df))) > 0) {
    names(df) <- sapply(names(df), function(x) strsplit(x, sep)[[1]][2])
  } 
  return(df)
}

softDown = function(GSE,file, overwrite=FALSE){
  if((file.exists(file) | file.exists(gsub('[.]gz', '', file))) & !overwrite){
    warning('this file already exists. not overwriting')
    return(FALSE)
  }
  download.file(paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
                       gsub('(((?<=GSE)([0-9]|[0-9][0-9]|[0-9][0-9][0-9]))|((?<=GSE.)[0-9][0-9][0-9])|((?<=GSE..)[0-9][0-9][0-9]))$|...$','nnn',GSE,perl = T),'/',
                       GSE,'/soft/',GSE,'_family.soft.gz'), file.path(file))
  return(TRUE)
}

ReadSoft <- function(fileName){ 
  con <- file(fileName,open="r")
  Meta <- c("Characteristics", "value")
  line <- readLines(con, n=1)
  while(length(line) > 0){
    line <- readLines(con, n=1)
    if(length(line) > 0) {
      if(grepl("(!Sample_characteristics|sample_id|Sample_title|!Sample_platform_id|right$|left$|Sample_source_name_ch1|!Sample_organism_ch1)", line, ignore.case = T)){
        temp <- gsub("!Sample_characteristics_ch1 = ?", "", line)
        temp <- strsplit(temp, "(:|=)")[[1]]
        line <- c(paste(temp[-length(temp)], collapse=":"), rev(temp)[1])
        Meta <- rbind(Meta, line)
      }
    }
  }
  close(con)
  Meta <- Meta[-1,]
  Meta <- apply(Meta, c(1,2), function(x) gsub("(^ | $|!)", "", x))
  Meta <- apply(Meta, c(1,2), function(x) gsub("(N[/ ]A)", NA, x))
  AllMetaVar <- unique(Meta[,1])
  MetaDF <- sapply(AllMetaVar, function(x) {
    for(char in paste0("\\",c("(", "[", "]", "^", "$", ".", "*", "+", "{", "}", "|"))){
      x <- gsub(char, paste0("\\\\",char), x )
    }
    Meta[grep(paste0("^",x,"$"), Meta[,1] ), 2]
  }, simplify = FALSE) %>%
    do.call(cbind, .) %>%
    as.data.frame(row.names = Meta[grep("sample_id", Meta[,1]),2])
  return(MetaDF)
}

DownloadCel <- function(celList, path){
  if(sum(grepl(paste0(path, "?"), list.dirs(full.names = FALSE))) == 0){
    dir.create(path)
  }
  celList2 <- celList[!celList %in% sapply(list.files(path), function(x) strsplit(x, "\\.")[[1]][1])]
  sapply(celList2, function(CelFile){
    print(paste0(CelFile, ": ", grep(CelFile, celList), "/", length(celList)))
    url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", CelFile)
    html <- paste(readLines(curl(url)), collapse="\n")
    GSElink <- grep('href=\"ftp', strsplit(html, " ")[[1]], value = T)
    GSElink <- GSElink[grepl("CEL", GSElink, ignore.case = TRUE)]
    GSElink <- strsplit(GSElink,'\\\"')[[1]][2]
    download.file(GSElink, destfile=paste0(path, "/", CelFile, ".CEL.gz"), method="curl", quiet=TRUE )
  })
}

DownloadExpFile <- function(study, path){
  if(sum(grepl(paste0(path, "?"), list.dirs(full.names = FALSE))) == 0){
    dir.create(path)
  }
  url <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", name)
  html <- paste(readLines(curl(url)), collapse="\n")
  GSElink <- grep('href=\"ftp', strsplit(html, " ")[[1]], value = T)
  GSElink <- GSElink[grepl("matrix", GSElink, ignore.case = TRUE)]
  GSElink <- strsplit(GSElink,'\\\"')[[1]][2]
  GSElink <- paste0(GSElink, name, "_series_matrix.txt.gz")
  download.file(GSElink, destfile=paste0(path, name, ".matrix.gz"), method = "curl", quiet=TRUE )
}

GetAnnoFiles <- function(platform){
  platform = as.character(platform)
  if(length(list.files(pattern = platform)) == 0){
    download.file(paste0("http://chibi.ubc.ca/microannots/", platform,  "_noParents.an.txt.gz"), destfile=paste0(platform, ".gz"))
  }
  if(length(list.files(pattern = platform)) == 1){
    warning("Using existing annotation file, consider updating")
  }
  Anno_file <- read.table(paste0(platform, ".gz"), comment="#", header=T, quote='"', sep="\t")
  return(Anno_file)
}
ScanData <- function(metaVar, scanVar){
  metaVar$ScanDate <- scanVar[match(metaVar$Filename, names(scanVar))]
  return(metaVar)
}

FindRep <- function(Metadata, char){
  names(Metadata) <- tolower(names(Metadata))
  char <- tolower(char)
  char <- char[char %in% names(Metadata)]
  AllChar <- Metadata %>% select_(.dots=char)
  Metadata$CharVec <- apply(AllChar, 1, function(x) paste0(x, collapse="_"))
  RepTable <- table(Metadata$CharVec)
  Replicates <- Metadata %>% filter(CharVec %in% names(RepTable[RepTable > 1])) %>% arrange(CharVec)
  Unique <- Metadata %>% filter(!CharVec %in% names(RepTable[RepTable > 1]))
  return(list(Replic = Replicates,
              Unique = Unique))
}

GetCommonName <- function(Metadata, Multi = Multi, char=char){
  if(tolower(Multi) == "yes"){
    SameSub <-  FindRep(Metadata=Metadata, char=char)
    UniqueName <- rbind(SameSub$Replic[!duplicated(SameSub$Replic$CharVec),],
                        SameSub$Unique) %>% arrange(profile)
    UniqueName$CommonName <- UniqueName$profile
    sampleSum <- table(UniqueName$CommonName)
    UniqueName$CommonName <- melt(sapply(names(sampleSum), function(x) {
      paste(x, seq(1:sum(UniqueName$CommonName == x)), sep="_")
      
    }))$value
    Metadata$CharVec <- sapply(Metadata$Filename,
                               function(x) rbind(SameSub$Replic, SameSub$Unique) %>%
                                 filter(filename == x) %>%
                                 .$CharVec) %>% unlist
    Metadata$CommonName <- UniqueName$CommonName[match(Metadata$CharVec, UniqueName$CharVec)]
    
  } else {
    
    Metadata$CommonName <- Metadata$Profile
    sampleSum <- table(Metadata$CommonName)
    
    Metadata$CommonName <- melt(sapply(names(sampleSum), function(x) {
      
      paste(x, seq(1:sum(Metadata$CommonName == x)), sep="_")
      
    }))$value
  }
  return(Metadata)
}

ModelAdj <- function(model, adj=data.frame(effect = "sex", adjValue=0)){
  type = class(model) %>% .[1]
  data <- model.frame(model)
  mod.matrix = as.matrix(model.matrix(model))
  for(i in 1:nrow(adj)){
    mod.matrix[,grep(adj$effect[i], colnames(mod.matrix), ignore.case=T)] <- adj$adjValue[i]
  }
  if(type=="lmerMod"){
    fixefVec = t(fixef(model)) %>% as.vector
    randAll <- ranef(model)
    randefList = sapply(names(randAll), function(randfcr) {
      sapply(data[[randfcr]], function(level){
        randAll[[randfcr]][level, "(Intercept)"]
      })
    }, simplify=FALSE)
    print(praise("${EXCLAMATION}! ${Adverb} ${created}"))
    randSum <- do.call(cbind, randefList) %>% rowSums
  } else {
    fixefVec = t(coef(model)) %>% as.vector
    randSum <- rep(0,length(resid))
  }
  fixSum = mod.matrix %*% fixefVec
  resid = resid(model)
  Final <- rowSums(cbind(fixSum, randSum, resid))
  return(Final)
}
