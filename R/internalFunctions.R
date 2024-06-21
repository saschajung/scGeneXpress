#' Data normalization UMI
#'
#' @param data_counts single cell matrix
#' @param log_transform log transformation
#' @param bg.model.pars background factors retrieved from scTransform (used when normalizing query)
#' @param data.type type of data: Ref (background) or query
#' @param file.id prefix to use when saving files
#' @param saveData save normalized data in rds object
#' @param return return normalized data from the function
#'
#' @return normalized data
scrb_norm<-function(data_counts,log_transform = T,cell.attr = NULL,bg.model.pars=NULL,data.type=c("Ref","query"),file.id="mydata",saveData=T,return=T){
  tictoc::tic("Data normalisation (UMI)")
  print("Normalising data")
  head(cell.attr)
  if(data.type == "Ref"){
    print("ref")
    data_norm<-DataNorm(data = data_counts,
                        log_transform = log_transform,
                        data.type = data.type,
                        file.id = file.id)
  }else if(data.type == "query"){
    print("query")
    #Keep genes found in query & background
    data_norm<-DataNorm(data = data_counts,
                        log_transform = log_transform,
                        data.type = data.type,
                        bg.model.pars = bg.model.pars,
                        cell.attr = cell.attr,
                        file.id = file.id)
  }
  tictoc::toc()
  print("saving data")
  if(saveData){
    if(data.type=="Ref"){
    dir_name=base::paste0("Ref_",file.id)
    }else if (data.type=="query"){
      dir_name=paste0("Query_",file.id)
    }
    base::dir.create(dir_name)
    base::saveRDS(data_norm,paste0(dir_name,"/",file.id,"_NormData.rds"))
  }
  if(return){
    return(data_norm)
  }
}

#' Data cleaning: removing genes not expressed, removing outliers & genes for which we cannot properly approximate distribution
#'
#' @param data_counts single cell matrix
#' @param dir.name directory where to save data
#' @param file.id prefix to use when saving data
#'
#' @return cleaned data
#' @author Céline Barlier, Sascha Jung
cleanData <- function(data_counts,dir.name,file.id,nmads = 3){
  print("Quality control")

  print("Removing cells with no counts")
  data_counts_cleaned <- data_counts[,Matrix::colSums(data_counts)>0]

  print("Removing low-quality cells")
  qcCells <- as.data.frame(scuttle::perCellQCMetrics(data_counts_cleaned))
  countsCells <- qcCells$sum
  base::names(countsCells) <- base::colnames(data_counts_cleaned)
  low.quality.cells <- scuttle::isOutlier(countsCells, nmads = nmads, log = TRUE)
  low.quality.cells <- base::names(low.quality.cells[base::which(low.quality.cells == TRUE)])
  data_counts_cleaned <- data_counts_cleaned[,!base::colnames(data_counts_cleaned) %in% low.quality.cells]

  print("Removing genes expressed in less than 10 cells")
  keep_genes <- Matrix::rowSums(data_counts_cleaned > 0) >= 10
  keep_genes <- base::names(base::which(keep_genes))
  data_counts_cleaned <- data_counts_cleaned[keep_genes,]

  output_list <- list("Counts" = data_counts_cleaned,
                      "QCMetrics" = qcCells)
  return(output_list)
}

#' Background sampling & threshold computation
#'
#' @param data_norm normalized single cell matrix
#' @param dir_name directory name where to save the data
#' @param file.id prefix to use when saving data
#' @param metadata metadata contains cells.id & cell.type
#' @param sample_n number of cells to use by group for sampling
#' @param parallelize run the code in parallel
#' @param ncores number of cores to use if parallelize = T
#' @param num.boot.sample number of bootstrap samples to run
#' @param fixseed seed used to generate ecdfBack
#'
#' @return background thresholds
#' @author Céline Barlier
GetBackground <- function(data_norm,dir_name,file.id,metadata,sample_n,parallelize,ncores,num.boot.sample,fixseed){
  #Variable init
  if(base::ncol(metadata) == 3){
    #Cell type level: two layers of information
    clusterids_vec = c("cell.type","tissue")
  }else{
    #Cell subtype or phenotypes level: one layer of information
    clusterids_vec = c("cell.type")
  }

  #Get original scaling factors and save them in file (used to scale queries later on) @TODO
  #Format data
  tictoc::tic("Preparing the background data")
  genes=base::rownames(data_norm)
  sf<-base::as.numeric(sparseMatrixStats::rowMaxs(data_norm))
  base::names(sf) <- genes
  #Save scaling factors
  utils::write.table(sf,file=base::paste0(dir_name,"/",file.id,"_sf.txt"),col.names = F)
  tictoc::toc()

  sf <- utils::read.table(base::paste0(dir_name,"/",file.id,"_sf.txt"),header = F)

  rnames <- base::rownames(data_norm)
  cnames <- base::colnames(data_norm)
  sData <- Matrix::summary(data_norm)
  sf$i <- base::match(sf$V1,base::rownames(data_norm))
  sData$x <- sData$x/sf$V2[base::match(sData$i,sf$i)]
  data_norm <- sparseMatrix(i=sData[,1], j=sData[,2], x=sData[,3],
                            dimnames = list(rnames,cnames))

  #Run bootstrapping of the: (sampling, scaling, thresholds computation) & get the resulting background thresholds
  back.thrs.ecdfs <- RunBootstrap(data_norm = data_norm,dir_name = dir_name,file.id = file.id,clusterids_vec = clusterids_vec,metadata = metadata,sample_n = sample_n,parallelize=parallelize,ncores=ncores,num.boot.sample = num.boot.sample, fixseed = fixseed)

  #Save data in folder if parameter saveData = True
  print("saving data")
  if(dir.exists(dir_name)==FALSE){
    dir.create(dir_name)
  }
  saveRDS(back.thrs.ecdfs,paste0(dir_name,"/",file.id,"_th_dist.rds"))

  return(back.thrs.ecdfs)
}

#' Identify genes with a biomdal distribution in the background
#'
#' @param data_norm normalized single cell matrix
#' @param dir_name directory name where to save the data
#' @param file.id prefix to use when saving data
#' @param metadata metadata contains cells.id & cell.type
#' @param sample_n number of cells to use by group for sampling
#' @param ncores number of cores to use if parallelize = Ts
#' @param parallelize run the code in parallel
#'
#' @return background thresholds
#' @author Céline Barlier
GetBimodalGenes <- function(runi,data_norm,dir_name,file.id,metadata,sample_n,ncores,parallelize=TRUE){

  if(parallelize){
    tictoc::tic("Preparing parallel running")
    cl = makePSOCKcluster(ncores, outfile="")
    registerDoParallel(cl)
    tictoc::toc()
  }

  #Variable init
  if(length(colnames(metadata)) == 3){
    #Cell type level: two layers of information
    clusterids_vec = c("cell.type","tissue")
  }else{
    #Cell subtype or phenotype level: one layer of information
    clusterids_vec = c("cell.type")
  }

  print("Identification of bimodal genes in the background data")

  print("Bootstrap: distributions compiling")
  tictoc::tic("Bootstrap: distributions compiling")
  bar=utils::txtProgressBar(min=0,max=100,style = 3)
  lg_list = foreach::foreach(n=1:100,.packages=c("foreach","Matrix","LaplacesDemon"),.export=c("ScaledData","sample.data","cells_exp","Background_sampling","RunBackComputations","SampleBackground","ScaleBackground"))%dopar%{
    options(future.rng.onMisuse="ignore")
    Sys.sleep(1)
    utils::setTxtProgressBar(bar,n)
    sampled.back <- SampleBackground(norm.data=data_norm,cell_clusterid_mat=metadata,ncells.sample=sample_n,clusterids_vec=clusterids_vec)

    #Scaling gene expression
    scaled.back <- ScaleBackground(sampled.back=sampled.back,file.id=file.id)
    processed.back <- scaled.back$data_processed

    #Return
    return(processed.back)
  }
  close(bar)
  tictoc::toc()

  #Close clusters & clean env
  if(parallelize){
    parallel::stopCluster(cl)
  }

  #Format distribution for each gene & identify bimodality
  print("Formatting of the distributions")
  ng <- length(lg_list[[1]])
  bar=utils::txtProgressBar(min=0,max=ng,style = 3)
  #for each gene
  lg_list_f <- lapply(X = 1:ng, function(n){
    utils::setTxtProgressBar(bar,n)
    #Bind all bootstrap runs
    gdist <- do.call("rbind",lapply(X = 1:length(lg_list),FUN = function(a){
      return(lg_list[[a]][[n]])
    }))
  })
  names(lg_list_f) = names(lg_list[[1]])
  close(bar)

  print("Identification of bimodal genes")
  bar=utils::txtProgressBar(min=0,max=length(lg_list_f),style = 3)
  lg <- sapply(seq(1,length(lg_list_f)),function(a){
    utils::setTxtProgressBar(bar,a)
    #Gene expression in the background
    gv <- lg_list_f[[a]]
    #Check if bimodal
    isbimo <- LaplacesDemon::is.bimodal(x=gv,min.size=0.1)
    #Return
    return(isbimo)
  })
  close(bar)
  #Add names
  names(lg) <- names(lg_list_f)
  #Remove non-bimodal genes (FALSE)
  lg <- lg[which(lg == TRUE)]

  #Save
  #saveRDS(names(lg),paste0(dir_name,"/",file.id,"_bimodGenes_",runi,".rds"))

  #Clean memory
  gc()

  return(names(lg))
}

#' Run Bootstrapping for the sampling, scalling & thresholds computation of the background
#'
#' @param data_norm normalized single cell matrix
#' @param file.id prefix to use when saving data
#' @param metadata metadata contains cells.id & cell.type
#' @param saveData save data in the folder
#' @param sample_n number of cells to use by group for sampling
#' @param parallelize run the code in parallel
#' @param ncores number of cores to use if parallelize = T
#' @param num.boot.sample number of bootstrap samples to run
#' @param fixseed seed used to generate ecdfBack
#'
#' @return bootstrapped thresholds & ecdfs functions
#' @author Céline Barlier
RunBootstrap <- function(data_norm,dir_name,file.id,clusterids_vec,metadata,sample_n,parallelize,ncores,num.boot.sample,fixseed){

  data_norm <- t(data_norm) #This is needed for the memory-efficient sampling approach!
  #Prepare cluster & core for the parallelization
  if(parallelize){
    tictoc::tic("Preparing parallel running")
    cl = makePSOCKcluster(ncores, outfile="")
    registerDoParallel(cl)
    tictoc::toc()
  }

  #Run the code with progress bar - get thresholds
  print("Generating thresholds distributions from the background data")

  tictoc::tic("Generating thresholds")
  bar=utils::txtProgressBar(min=0,max=num.boot.sample,style = 3)
  th_list = foreach::foreach(n=1:num.boot.sample,.packages=c("foreach","Matrix"),.export=c(".optLowFun",".optHighFun","swap_th","minimizeRectangle","createThresholdDist","ScaledData","sample.data","cells_exp","Background_sampling","RunBackComputations","SampleBackground","ScaleBackground","GetThresholds"))%dopar%{
    options(future.rng.onMisuse="ignore")
    Sys.sleep(1)
    utils::setTxtProgressBar(bar,n)
    return(RunBackComputations(norm.data=data_norm,cell_clusterid_mat=metadata,ncells.sample=sample_n,clusterids_vec=clusterids_vec,file.id=file.id))
  }
  close(bar)
  tictoc::toc()

  #Combine lists to have thresholds by gene
  print("Formatting background thresholds")
  tictoc::tic("Formatting background thresholds")
  ng <- length(th_list[[1]])
  bar=utils::txtProgressBar(min=0,max=ng,style = 3)
  #for each gene
  thrs_list <- lapply(X = 1:ng, function(n){
    utils::setTxtProgressBar(bar,n)
    #Bind all bootstrap runs
    gthrs <- do.call("rbind",lapply(X = 1:length(th_list),FUN = function(a){
      return(th_list[[a]][[n]])
    }))
  })
  names(thrs_list) = names(th_list[[1]])
  close(bar)
  tictoc::toc()

  #Close clusters & clean env
  parallel::stopCluster(cl)

  return(thrs_list)
}

#' Function called by RunBootstrap() to sample, scale & compute the background thresholds
#'
#' @param norm.data normalized single cell matrix
#' @param cell_clusterid_mat metadata
#' @param ncells.sample number of cells to use by group for the sampling
#' @param clusterids_vec metadata column name with the cell type info
#' @param file.id prefix for the file to save
#'
#' @return sampled background
#' @author Modified by Céline Barlier
RunBackComputations <- function(norm.data,cell_clusterid_mat,sf,ncells.sample,clusterids_vec,file.id){

  #Sampling cells
  sampled.back <- SampleBackground(norm.data=norm.data,cell_clusterid_mat=cell_clusterid_mat,ncells.sample=ncells.sample,clusterids_vec=clusterids_vec,file.id=file.id)

  #Scaling gene expression
  #scaled.back <- ScaleBackground(sampled.back=sampled.back,file.id=file.id)
  #processed.back <- scaled.back$data_processed

  #Get thresholds
  thrs.back <- GetThresholds(processed.back=sampled.back,file.id=file.id)

  return(thrs.back)
}

#' Sampling of the background
#'
#' @param norm.data normalized single cell matrix
#' @param cell_clusterid_mat metadata
#' @param ncells.sample number of cells to use by group for the sampling
#' @param clusterids_vec metadata column name with the cell type info
#' @param file.id prefix for the file to save
#'
#' @return sampled background
#' @author Sascha Jung, Céline Barlier
SampleBackground <- function(norm.data,cell_clusterid_mat,ncells.sample,clusterids_vec,file.id){

  #In contrast to the previous version (below), this assumes that the data is already transposed!
  norm.data.list <- list()
  norm.data.list <- lapply(2:length(norm.data@p),function(i){
    sample(norm.data@x[(norm.data@p[i-1]+1):norm.data@p[i]],ncells.sample,replace = T)
  })
  names(norm.data.list) <- colnames(norm.data)
	
  #norm.data <- t(norm.data)
  #m2<-summary(norm.data)
  #norm.data.list<-split(m2$x,colnames(norm.data)[m2$j])

  #norm.data.list <- lapply(norm.data.list,function(x){
  #  sample(x,ncells.sample,replace = T)
  #})

  return(norm.data.list)
}

#' Scaling of the background
#'
#' @param sampled.back sampled single cell matrix
#' @param file.id prefix for the file to save
#'
#' @return sampled background
#' @author Modified by Céline Barlier
ScaleBackground <- function(sampled.back,file.id){
  scaled.back <- ScaledData(sampled.back)
  return(scaled.back)
}

#' Get thresholds of the background
#'
#' @param processed.back sampled & sclaed single cell matrix
#' @param file.id prefix for the file to save
#'
#' @return sampled background
#' @author Modified by Céline Barlier
GetThresholds <- function(processed.back,file.id){
  check=names(which(sapply(processed.back,length)<1))
  if(length(check)>0){
    processed.back=processed.back[setdiff(names(processed.back),check)]
  }

  th_list = lapply(processed.back, function(pbx){
    createThresholdDist(pbx)
  })
  names(th_list)=names(processed.back)

  return(th_list)
}

#' Get ecdfs of the background
#'
#' @param processed.back sampled & sclaed single cell matrix
#' @param file.id prefix for the file to save
#'
#' @return sampled background
#' @author Modified by Céline Barlier
GetEcdfs <- function(processed.back,file.id){
  check=names(which(sapply(processed.back,length)<1))
  if(length(check)>0){
    processed.back=processed.back[setdiff(names(processed.back),check)]
  }

  ecdf_list = lapply(processed.back, function(pbx){
    createEcdfsBack(pbx)
  })
  names(ecdf_list)=names(processed.back)

  return(ecdf_list)
}

#' Prepare query data: clean, normalize & scale
#'
#' @param query_dat query data
#' @param Ref_data background reference
#' @param dir directory used to save the results
#' @param sfback scaling factors
#' @param file.id prefix to use when saving files
#' @param ncores number of cores to use to run the code
#' @param bg.model.pars background factors retrieved from scTransform (used when normalizing query)
#' @param normalise normalize data
#' @param scale_data scale data
#' @param save data in rds object
#' @param return.data return results
#'
#' @return prepared query data
PrepQueryData<-function(query_dat,
                        Ref_data,
                        dir,
                        sfback=NA,
                        file.id="query",
                        ncores,
                        bg.model.pars,
                        normalise=T,
                        scale_data=T,
                        saveData=T,
                        return.data=T){
  scaling.factors.file = list.files(Ref_data,pattern = "sf",full.names = T)
  dir_name=dir
  if(dir.exists(dir)==FALSE){
    dir.create(dir)
  }
  log.file=paste0(dir_name,"/",file.id,"_log.txt")
  cat("Log file",file=log.file,sep="\n",append=T)
  cat(paste("Find data in",paste0(getwd(),"/",dir_name)),file=log.file,sep="\n",append = T)
  cat("Preparing data: Normalisation and Scaling",file=log.file,sep="\n",append = T)
  tic.clearlog()
  print("Cleaning & QC of the data")
  query_dat <- cleanData(data_counts=query_dat,
                         dir.name=paste0(getwd(),"/",dir_name),
                         file.id=file.id)

  query_dat$QCMetrics$log_umi <- log10(query_dat$QCMetrics$sum)
  query_dat$QCMetrics$detected <- NULL
  query_dat$QCMetrics <- query_dat$QCMetrics[,c("sum","log_umi")]
  colnames(query_dat$QCMetrics) <- c("umi","log_umi")

  ##normalisation
  if(normalise){
	  print("Normalising data (UMI): ScTransform")
	  tictoc::tic("Data Normalisation")
	  query_Normdata=scrb_norm(data_counts = query_dat$Counts,
	                           file.id = file.id,
	                           data.type="query",
	                           saveData = T,
	                           return = T,
	                           bg.model.pars = bg.model.pars,
	                           cell.attr = query_dat$QCMetrics)
	    tictoc::toc(log=T)
  }else{
    query_Normdata=query_dat$Counts
  }

  ##scaling
  #If no pre-compiled background
  if(class(sfback) != "data.frame"){
    scaling.factors.file = list.files(Ref_data,pattern = "sf",full.names = T)
    sf=read.delim(scaling.factors.file,sep=" ",row.names = 1,header = F)
  }else{
    sf=sfback
  }
  if(scale_data){
    genes_common<-intersect(rownames(query_Normdata),rownames(sf))
    print(paste("Number of genes in the query with corresponding reference distributions: ",length(genes_common)))
    print("Scaling query data with scaling factors obtained from the background")
    tictoc::tic("Scaling data")
    query_scaled_data=scaleMaxGexp(query_Normdata,sf,genes_common)
    tictoc::toc(log=T)
    if(saveData){
      print("saving scaled data")
      saveRDS(query_scaled_data,file=paste0(dir_name,"/",file.id,"_ScaledData.rds"))
    }
  }
  log=tic.log(format = T)
  tic.clearlog()
  cat("Time log",file=log.file,sep="\n",append = T)
  cat(unlist(log),file=log.file,sep="\n",append = T)
  cat("End Time log",file=log.file,sep="\n",append = T)
  if(return.data){
  return(query_scaled_data)
  }
}

#' Discretize gene expression in three levels of expression: 0 (not expressed), 0.5 (medium level), 1 (high level of expression)
#'
#' @param query_ScaledData prepared query data
#' @param Ref_data background reference
#' @param thsback background thresholds
#' @param th threshold cutoffs
#' @param th_i threshold cutoffs
#' @param file.id prefix to use when saving files
#' @param return.data return results
#' @param save data in rds object
#'
#' @return discretized gene expression matrix
Discretise_data<-function(query_ScaledData,Ref_data,thsback=NA,th=0.05,th_i=0.05,file.id="myquery",return.data=T,saveData=T){
  dir_name=paste0("Query_",file.id)
  if(dir.exists(dir_name)==FALSE){
    dir.create(dir_name)
  }

  log.file=paste0(dir_name,"/",file.id,"_log.txt")

  cat("Discretising data",file=log.file,sep="\n",append = T)
  tic.clearlog()

  if(length(thsback)==1){
    threshold.dist.file=list.files(Ref_data,pattern = "_dist",full.names = T)
    threshold.dist=readRDS(threshold.dist.file)
  }else{
    threshold.dist = thsback
  }

  ##computing pvals
  genes_common<-intersect(rownames(query_ScaledData),names(threshold.dist))
  tictoc::tic("computed q values")
  query_pvals=compute_pval(th_list = threshold.dist,gexp_mat = query_ScaledData,genes = genes_common,p_correct = F)
  tictoc::toc(log=T)
  saveRDS(object = query_pvals,file = paste0(dir_name,"/",file.id,"_pvals.rds"))

  ##discretising gene expression
  tictoc::tic("Discretised data")
  query_discretised = classify(query_pvals,th,th_i)
  tictoc::toc(log=T)
  data=list(discretised_data=query_discretised,significance=query_pvals)
  if(saveData){
    print("saving discretised data")
    saveRDS(data,file=paste0(dir_name,"/",file.id,"_DiscretisedData.rds"))
  }
  log=tic.log(format = T)
  tic.clearlog()
  cat("Time log",file=log.file,sep="\n",append = T)
  cat(unlist(log),file=log.file,sep="\n",append = T)
  cat("End Time log",file=log.file,sep="\n",append = T)
  if(return.data){
  return(data)
  }
}

#' Function to quantify gene level at the population level and identify key identity genes
#'
#' @param discrete.mat discretized matrix
#'
#' @return data.frame of gene levels frequency
#@Modified by Céline Barlier
Create.sig.frame <- function(discrete.mat){

  discrete.mat[is.nan(discrete.mat)] <- NA

  p <- sum(discrete.mat > 0,na.rm = TRUE)/c(nrow(discrete.mat) * ncol(discrete.mat))
  n <- ncol(discrete.mat)

  mean = n * p
  sd = sqrt(n * p * (1-p))

  z2 = (sd * 2) + mean

  sig.frame <- do.call("rbind",lapply(X = row.names(discrete.mat),FUN = function(g){

    g.discrete.mat <- discrete.mat[g,]

    if(all(is.na(g.discrete.mat))){
      cat(" All values NA\n")
      return(NULL)
    }else{

      n.0 <- length(g.discrete.mat[g.discrete.mat == 0 & !is.na(g.discrete.mat)])
      n.0.5 <- length(g.discrete.mat[g.discrete.mat == 0.5 & !is.na(g.discrete.mat)])
      n.1 <- length(g.discrete.mat[g.discrete.mat == 1 & !is.na(g.discrete.mat)])
      not.0 <- n.0.5+n.1
      not.0.perc <- (not.0)/(n.0+n.0.5+n.1)
      row <- cbind.data.frame(g,n.0,n.0.5,n.1,not.0,not.0.perc,stringsAsFactors = F)
      colnames(row)[1] <- "gene"
      return(row)
    }
  }))
  row.names(sig.frame) <- sig.frame$gene
  #sig.frame <- sig.frame[which(sig.frame$not.0 > z2),]
  sig.frame <- sig.frame[order(sig.frame$not.0.perc,decreasing = T),]
  return(list(sig.frame,z2))
}

#' Function to determine the gene expression level of key identity genes
#'
#' @param tbl sig.frame data.frame obtained from Create.sig.frame()
#' @param zs z-score computed in Create.sig.frame()
#'
#' @return key identity genes data.frame (with gene quantification level)
#' @author Céline Barlier
discrete.expr.levels <- function(tbl,zs){

  tbl$expr.level <- do.call("c",lapply(X = 1:nrow(tbl),function(i){

    #Zscore filtering
    if(tbl[i,"not.0"] <= zs){
      #Significantly lowly expressed
      return("low")
    }else{
      #Significantly mediumly or highly expressed: take the maximum value to determine
      n.0.5 <- as.numeric(tbl[i,"n.0.5"])
      n.1 <- as.numeric(tbl[i,"n.1"])
      #If got the exact same value we cannot determine > undefined
      if(n.0.5 == n.1){
        return("not.significant")
      }else{
        #If n0.5 > n.1 = medium
        if(n.0.5 > n.1){
          return("medium")
        }else{
          #n.1 > n0.5 = high
          return("high")
        }
      }
    }
  }))
  return(tbl)
}

#' Function to get identity genes: high level of expression (from discrete.expr.levels)
#'
#' @param tbl sig.frame data.frame obtained from Create.sig.frame()
#' @param bimodg genes with a bimodal distribution in the background (used to identify "unique" mediumly expressed genes)
#'
#' @return data.frame of identity genes and their level (high/medium)
#' @author Céline Barlier
get.identity.genes <- function(tbl,bimodg){

  #1. Add the high level of expression genes ("unique/specific" ones according to the background)
  kigHigh <- tbl$gene[which(tbl$expr.level == "high")]

  #Data.frame
  dfkig <- data.frame("IdentityGene"=kigHigh,"Level"=rep("high",length(kigHigh)))

  #2. Add the "unique/specific" medium genes (= the gene is predicted as medium in the population & its distribution in the background is bimodale)
  kigMed <- tbl$gene[which(tbl$expr.level == "medium")]
  kigMed <- kigMed[which(kigMed %in% bimodg)]

  #Add to data.frame if any
  if(length(kigMed)>0){
    dfkig <- rbind(dfkig,data.frame("IdentityGene"=kigMed,"Level"=rep("medium",length(kigMed))))
  }

  #Return
	return(dfkig)
}

#' Data Normalisation
#'
#' @param data single cell data to normalize
#' @param log_transform log transformation
#' @param data.type type of data (Ref or query)
#' @param bg.model.pars background factors retrieved from scTransform (used when normalizing query)
#' @param file.id presfix to use when saving files
#'
#' @return normalized data
DataNorm<-function(data,log_transform=c(TRUE,FALSE),cell.attr,scale_factor = NA,data.type,bg.model.pars=NULL,file.id=NULL){
  if(data.type == "Ref"){
    if(file.exists(paste0("./Ref_",file.id,"/",file.id,"_bg_model_pars.rds"))){
      vst.out <- readRDS(paste0("./Ref_",file.id,"/",file.id,"_bg_model_pars.rds"))
    }else{
    vst.out <-sctransform::vst(umi = data,
                               return_corrected_umi = FALSE,
                               n_genes=NULL,
                               n_cells = base::ncol(data),
                               min_cells=5,
                               return_cell_attr = TRUE,
                               return_gene_attr = TRUE,
                               scale_factor = scale_factor,
                               vst.flavor = "v2",
                               residual_type = "none")
    saveRDS(object = vst.out,file = paste0("./Ref_",file.id,"/",file.id,"_bg_model_pars.rds"))
    }
    cell.attr <- vst.out$cell_attr
  }else if(data.type == "query"){
    vst.out <- bg.model.pars
    data <- data[which(rownames(data) %in% rownames(vst.out$model_pars)),]
    cell.attr <- cell.attr[which(rownames(cell.attr) %in% colnames(data)),]
  }
  dat_norm <- sctransform::correct_counts(x = vst.out,
                                          umi = data,
                                          cell_attr = cell.attr,
                                          scale_factor = vst.out$arguments$scale_factor)

  if(log_transform){
    return(base::log1p(dat_norm))
  }else{
    return(dat_norm)
  }
}

#' Sample cells
#'
#' @param cells_dat vector of cell names/UMI
#' @param cellid named vector (names=cell name/UMI, elements=clusterid-combined id)
#' @param ncells.sample number of cells to sample
#' @param nlevel number of cluster id levels (e.g Tissue , cell type, cell subtype)
#'
#' @return sampled cells
sample.data<-function(cells_dat,cellid,ncells.sample,nlevels){

  #If only one level provided (cell type or cell subtype)
  if(nlevels==1){
    #Init var
    set_samp=cellid
    cell.id.split=set_samp

    #Split cells by cell type
    set=split(set_samp,as.factor(cell.id.split))

    #For each cell type
    set_cells=unlist(sapply(set,function(x){
      pick=intersect(names(cells_dat),names(x))
      if(length(pick)>=1&&length(x)>=50){
        sub_cell= sample(x = pick,size = ncells.sample,replace=T)
      }else{
        sub_cell=NA
      }
      return(sub_cell)
    }))
  }else{
    #2 levels provided: tissue & cell type
    set_samp=cellid
    cell.id.split=sapply(set_samp,function(x){paste(strsplit(x,"[.]")[[1]][1],collapse = ".")}) #CT

    #Split cells by cell type
    set=split(set_samp,as.factor(cell.id.split))

    #For each cell type
    set_cells=unlist(sapply(set,function(x){

      #Split by tissue
      cell.id.split.tissue=sapply(x,function(y){paste(strsplit(y,"[.]")[[1]][2],collapse = ".")})
      set.tissue=split(x,as.factor(cell.id.split.tissue))
      #Number of cells by tissue
      nbc = unlist(lapply(set.tissue,function(b){length(b)}))
      #How many to pick up by tissue - equally distributed
      n=round(ncells.sample/length(set.tissue))
      #Keep the ones >= n
      set.tissue=set.tissue[which(as.numeric(nbc) >= n)]
      #Size to sample
      sn=round(ncells.sample/length(set.tissue))

      #For each tissue
      set_cells_tis=unlist(sapply(set.tissue,function(a){
        picktis=intersect(names(cells_dat),names(a))
        if(length(picktis)>=1&&length(a)>=sn){
          sub_cell=sample(x = picktis,size = sn,replace=T)
        }else{
          sub_cell=NA
        }
        return(sub_cell)
      }))
      set_cells_tis=set_cells_tis[which(is.na(set_cells_tis)==FALSE)]
      #If at least 50 cells sampled
      if(length(set_cells_tis)>=50){
        #If less than 100 cells : oversample, if more than 100: downsample
        if(length(set_cells_tis) < 100 || length(set_cells_tis) > 100){
          set_cells_tis=sample(x=set_cells_tis,size=ncells.sample,replace=T)
        }
      }else{
        set_cells_tis=NA
      }
      return(set_cells_tis)
    }))
  }
  if(is.null(length(set_cells))){
    return(NULL)
  }else{
    set_cells=set_cells[which(is.na(set_cells)==FALSE)]
    set_samp=set_samp[set_cells]
    x=cells_dat[names(set_samp)]
    names(x)=names(set_samp)
    return(x)
  }
}

#' Filter cells with non zero gene expression
#'
#' @param data_matrix matrix genes x UMIs
#' @param min.cells.gexp minimum of cells genes need to be expressed in
#'
#' @return list of genes with names of cells (UMIs) that have non zero gene exp values
cells_exp<-function(data_matrix,min.cells.gexp){
  #In case matrix is too big, split by 10k rows to avoid chmold errors
  dimmtx <- dim(data_matrix)
  #If more than 10k genes & more than 50k cells
  if(dimmtx[1]>10000 & dimmtx[2]>50000){
    #Compute number of chunks: every 5k genes
    x <- seq(1,nrow(data_matrix))
    n <- round(nrow(data_matrix)/5000)
    chunks <- split(x, sort(x%%n))
    #Init list
    data_nonzero_list=apply(data_matrix[chunks[[1]],],1,function(x){
        return(x[x>0])
    })
    #Complete list
    for(i in seq(2,length(chunks))){
      runlist=apply(data_matrix[chunks[[i]],],1,function(x){
        return(x[x>0])
      })
      #Bind lists
      data_nonzero_list = c(data_nonzero_list,runlist)
      rm(runlist)
    }
  }else{
    data_nonzero_list=apply(data_matrix,1,function(x){
      return(x[x>0])
    })
  }

  len=sapply(data_nonzero_list,length)
  return(data_nonzero_list[len>min.cells.gexp])
}

#' Background sampling
#'
#' @param counts_norm_mat normalised umi counts matrix-- genes X UMIs
#' @param cell_clusterid_mat matrix with columns --"cell.name(UMI)"--same as colnames of counts matrix,"cluster.id","Tissue"
#' @param ncells.sample number of cells to sample
#' @param clusterids_vec vector with elements corresponding column names in cell_clusterid_mat with the correct heirarchiel order of clusterids (eg. Tissue, cell type, cell subtype)
#' @param min.cells.gexp minimum of cells genes need to be expressed in
#'
#' @return sampled background
Background_sampling<-function(counts_norm_mat,cell_clusterid_mat,ncells.sample=100,clusterids_vec="cell.type",min.cells.gexp=10){

  nlevels=length(clusterids_vec)
  cellid<-apply(cell_clusterid_mat,1,function(x){paste(x[clusterids_vec],collapse = ".")})
  names(cellid)=cell_clusterid_mat$cell.id
  cells_sample=cells_exp(counts_norm_mat,min.cells.gexp=min.cells.gexp)
  genes=names(cells_sample)

  ids_samples_byGene=lapply(cells_sample, function(csl){
    return(sample.data(cells_dat = csl,cellid = cellid,ncells.sample =  ncells.sample, nlevels = nlevels))
  })
  names(ids_samples_byGene)=genes
  #Remove NULL (genes that did not pass the criteria)
  ids_samples_byGene=ids_samples_byGene[!sapply(ids_samples_byGene,is.null)]

  return(ids_samples_byGene)
}

#' Scaling of data
#'
#' @param data single cell matrix to scale with maximum gene expression
#'
#' @return scaled matrix
ScaledData<-function(data){
  #Scaling factors
  sf<-sapply(data,max)
  #Format & scale data
  genes=names(data)
  data_new<-lapply(1:length(sf),function(x){
    return(data[[x]]/sf[x])
  })
  #Return
  names(data_new)=genes
  return(list(data_processed=data_new,scaling_factors=sf))
}

#' Called by ScaledData()
#'
#' @param data_mat single cell matrix
#' @param scalingfactors scaling factors
#' @param genes gene names
#'
#' @return scaled expression
scaleMaxGexp<-function(data_mat,scalingfactors,genes){
  #bar=utils::txtProgressBar(min=0,max=length(genes),style = 3)
  data_s=foreach(n=1:length(genes))%do%{
    #utils::setTxtProgressBar(bar,n)
    data=data_mat[genes[n],]/scalingfactors[genes[n],]
    return(data)
  }
  #close(bar)
  names(data_s)=genes
  data_s=do.call(rbind,data_s)
  return(data_s)
}

#' Optimization function to compute lower-bound thresholds
.optLowFun <- function(x,empcdf,maxVal){
  return(abs(x* (1-empcdf(x))))
}

#' Optimization function to compute upper-bound thresholds
.optHighFun <- function(x,empcdf,maxVal){
  return(abs(empcdf(x) * (maxVal-x)))
}

#' Function to compute thresholds
minimizeRectangle <- function(empCDF,maxVal,minVal){
  lowThr <- optimize(f = .optLowFun, interval = c(minVal,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum

  highThr<- optimize(f = .optHighFun, interval = c(minVal,maxVal),  empcdf = empCDF, maxVal = maxVal, maximum = TRUE)$maximum

  return(c(lowThr,highThr))
}

#' Verify consistency of computed thresholds
#' Check that lower gexp value are lower thresholds & higher gexp value are upper thresholds
#'
#' @param th_dist thresholds
#'
#' @return thresolds
swap_th<-function(th_dist){
  mat=apply(th_dist,1,function(x){
    if(x[1]>x[2]){
      temp=x[1]
      x[1]=x[2]
      x[2]=temp
      return(c(x[1],x[2]))
    }else{
      return(c(x[1],x[2]))
    }
  })
  return(t(mat))
}

#' Generate threshold distributions
#'
#' @param gexp gene expression vector
#' @param swap parameter to ensure consistency of lower & upper-bounds
#'
#' @return thresholds
#Modified by Céline Barlier
createThresholdDist <- function(gexp,swap=TRUE){
  gexp=as.numeric(gexp)

  if(length(unique(gexp))==1){
    thrs<-c(0,unique(gexp))
  }else{
    thrs<-minimizeRectangle(ecdf(gexp),max(gexp,na.rm=T),min(gexp,na.rm=T))
  }
  thrsdf <- data.frame("Lower"=thrs[1],"Upper"=thrs[2])
  #Define lower gexp value to be lower-threshold and higher gexp to be upper-threshold
  if(swap){
    thrsdf=swap_th(thrsdf)
  }
  return(thrsdf)
}

#' Generate ecdf function for each gene
#'
#' @param gexp gene expression vector
#'
#' @return ecdfs
#' @author Céline Barlier
createEcdfsBack <- function(gexp,swap=TRUE){
  gexp=as.numeric(gexp)
  return(ecdf(gexp))
}

#' Compute p-values
#'
#' @param th_list list of background thresholds
#' @param gexp_mat gene expression matrix (query)
#' @param genes genes for which p-values will be computed
#' @param p_correct perform multiple correction
#'
#' @return p-values (corrected if p_correct = T)
compute_pval<-function(th_list,gexp_mat,genes,p_correct=T){
  print("computing p values from threshold distributions")
  bar=utils::txtProgressBar(min=0,max=length(genes),style = 3)
  p_vals=foreach(n=1:length(genes))%do%{
    utils::setTxtProgressBar(bar,n)
    x=genes[n]
    p_vals_exp=1-(ecdf(th_list[[x]][,"Upper"])(gexp_mat[x,]))
    p_vals_nexp=ecdf(th_list[[x]][,"Lower"])(gexp_mat[x,])
    return(list(p_vals_exp=p_vals_exp,p_vals_nexp=p_vals_nexp))
  }
  close(bar)
  if(p_correct){
    print("Correcting for multiple testing: computing q values")
    sig_vals=lapply(p_vals,function(x){
      q_vals_exp=stats::p.adjust(x$p_vals_exp,method = "BH")
      q_vals_nexp=stats::p.adjust(x$p_vals_nexp,method = "BH")
      return(list(sig_exp=q_vals_exp,sig_nexp=q_vals_nexp))
    })
  }else{
    sig_vals=lapply(p_vals,function(x){
      return(list(sig_exp=x$p_vals_exp,sig_nexp=x$p_vals_nexp))
  })
  }
  names(sig_vals)=genes
  sig_exp=do.call("rbind",lapply(sig_vals,function(x){x$sig_exp}))
  sig_nexp=do.call("rbind",lapply(sig_vals,function(x){x$sig_nexp}))
  colnames(sig_exp)=colnames(gexp_mat)
  colnames(sig_nexp)=colnames(gexp_mat)
  return(list(sig_exp=sig_exp,sig_nexp=sig_nexp))
}

#' Classify gene activity
classify<-function(sig_vals,th,th_i){
  print("Discretising gene expression values")
  discrete_exp<-matrix(NaN,nrow=dim(sig_vals[[1]])[1],ncol=dim(sig_vals[[1]])[2])
  rownames(discrete_exp)=rownames(sig_vals[[1]])
  colnames(discrete_exp)=colnames(sig_vals[[1]])
  bar=utils::txtProgressBar(min=0,max=dim(sig_vals[[1]])[1],style = 3)
  for(i in 1:dim(discrete_exp)[1]){
    utils::setTxtProgressBar(bar,i)
    for(j in 1:dim(discrete_exp)[2]){
      if(is.nan(sig_vals[["sig_exp"]][i,j])){
        discrete_exp[i,j]=0
      }
      if(is.nan(sig_vals[["sig_nexp"]][i,j])){
        discrete_exp[i,j]=0
      }else if(sig_vals[["sig_exp"]][i,j]<=th){
        discrete_exp[i,j]=1
      }else if(sig_vals[["sig_nexp"]][i,j]<=th){
        discrete_exp[i,j]=0
      }else if(sig_vals[["sig_nexp"]][i,j]>1-th_i&&sig_vals[["sig_exp"]][i,j]>th_i){
        discrete_exp[i,j]=0.5
      }
    }
  }
  close(bar)
  return(discrete_exp)
}

