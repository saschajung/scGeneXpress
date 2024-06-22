#' scGeneXpress is a computational method for detecting expressed genes in scRNAseq samples
#'
#' @param data single cell expression matrix
#' @param metadata corresponding metadata (two columns: cell names, cell type)
#' @param org organism (human or mouse)
#' @param dir directory used to save the results
#' @param file.name prefix for output files and folder
#' @param precBack precompiled background path to ScRbRef
#' @param discretize gene quantification in three levels of cell types
#' @param sig.frames construction of significant genes data frame
#' @param ncores number of cores to use to run the code
#' @param fixseed seed used to discretize cell types (reproducbility)
#'
#' @import foreach
#' @import parallel
#' @import doParallel
#' @import reshape2
#' @import tictoc
#' @import SingleCellExperiment
#' @import Matrix
#' @import scater
#' @import stringr
#' @import qlcMatrix
#' @import sctransform
## @import org.Mm.eg.db
## @import org.Hs.eg.db
#' @import scuttle
#' @import LaplacesDemon
#'
#' @export
#'
#' @return NULL, all results are accessible in the dir specified
#' @author Céline Barlier, Sascha Jung
run_scGeneXpress <- function(data,metadata,org="human",dir,file.name,precBack=NA,discretize=TRUE,sig.frames=TRUE,ncores=detectCores()-2,fixseed=1234){

  if(base::endsWith(dir,"/")){
    dir <- base::substr(dir, 1, base::nchar(dir)-1)
  }

  # Changing column names of metadata for consistency
  if(base::ncol(metadata) == 3){
    base::colnames(metadata) <- c("cell.id","tissue","cell.type")
  }else{
    base::colnames(metadata) <- c("cell.id","cell.type")
  }
  base::rownames(metadata) <- metadata$cell.id

  # Transforming matrix to sparse matrix if needed
  if(base::is.matrix(data) | base::is.data.frame(data)){
    data <- as(as.matrix(data), "dgCMatrix")
  }

  # Identifying common cells between gene expression data frame and metadata - cleaning data and metadata
  common.cells <- intersect(metadata$cell.id,base::colnames(data))
  metadata <- metadata[base::which(metadata$cell.id %in% common.cells),]
  data <- data[,base::colnames(data) %in% common.cells]

  ################### BACKGROUND CONSTRUCTION/SELECTION #######################

  #If not pre-compiled background specified, build background
  if(is.na(precBack)){

    # Calculating number of cells of different cell types in data
    cell.type.freq <- as.data.frame(base::table(metadata$cell.type))
    # Need at least 50 cells by cell type
    cell.type.to.keep <- base::as.character(cell.type.freq[which(cell.type.freq$Freq >= 50),"Var1"])
    metadataBack <- metadata[base::which(metadata$cell.type %in% cell.type.to.keep),]
    dataBack <- data[,base::as.character(metadataBack$cell.id)]

    if(base::length(cell.type.to.keep) < 2){

      stop("Data contains less than 2 celltypes with at least 50 cells. Background reference cannot be created. Aborting...")

    }else{

      removed.cell.type <- base::setdiff(cell.type.freq$Var1,cell.type.to.keep)
      if(base::length(removed.cell.type) == 0){
        removed.cell.type <- "None"
      }

      #Create directory for results
      if(!base::dir.exists(paths = base::paste0(dir,"/Ref_",file.name))){
        base::dir.create(base::paste0(dir,"/Ref_",file.name))
      }

      #Create log file
      log.file.name <- base::paste0(dir,"/Ref_",file.name,"/Log_file_",base::paste(base::unlist(base::strsplit(x = base::as.character(base::Sys.time()),split = " ",fixed = 2)),collapse = "_"),".txt")
      write(x = base::paste0("Cell types used for Reference background: ",base::paste(cell.type.to.keep,collapse = ",")),file = log.file.name)
      write(x = base::paste0("Cell types removed: ",base::paste(removed.cell.type,collapse = ",")),file = log.file.name,append = T)

      setwd(dir)
      #Cleaning & Normalizing background data
      if(!base::file.exists(base::paste0("./Ref_",file.name,"/",file.name,"_NormData.rds"))){
        clean.data <- base::suppressWarnings(cleanData(data_counts = dataBack,dir.name=base::paste(dir,paste("/Ref_",file.name,sep=""),sep=""),file.id = file.name))
        norm.data <- base::suppressWarnings(scrb_norm(data_counts = clean.data$Counts,file.id = file.name,data.type = "Ref",saveData = T,return = T))
        rm(clean.data)
        invisible(gc())
      }else{
        norm.data <- base::readRDS(base::paste0("./Ref_",file.name,"/",file.name,"_NormData.rds"))
      }

      #Sampling, scaling & getting background thresholds
      if(!file.exists(paste0("./Ref_",file.name,"/",file.name,"_th_dist.rds"))){
        thrs.data <- GetBackground(data_norm=norm.data,dir_name=base::paste0(dir,"/Ref_",file.name),file.id=file.name,metadata=metadataBack,sample_n = 100,parallelize=T,ncores=ncores,num.boot.sample = 1000,fixseed = fixseed)
      }else{
        thrs.data <- readRDS(paste0("./Ref_",file.name,"/",file.name,"_th_dist.rds"))
      }

      #Identify biomdal genes in the background
      # if(!file.exists(paste0("./Ref_",file.name,"/",file.name,"_bimod.rds"))){
      #   #3 runs to ensure stability
      #   bgl <- list()
      #   for(i in seq(1,3)){
      #     print(paste0("Run ",i,"/3"))
      #     bgl[[i]] <- GetBimodalGenes(runi=i,data_norm=norm.data,dir_name=paste0(dir,"/Ref_",file.name),file.id=file.name,metadata=metadataBack,sample_n = 100,ncores=ncores)
      #   }
      # 	bimod.genes <- Reduce(intersect,bgl)
      #   saveRDS(bimod.genes,paste0("./Ref_",file.name,"/",file.name,"_bimod.rds"))
      # }else{
      # 	bimod.genes <- readRDS(paste0("./Ref_",file.name,"/",file.name,"_bimod.rds"))
      # }

      setwd(paste0("./Ref_",file.name))
      #Normalization parameters used for the background (the algorithm will use the same to normalize the queries)
      bg.model.pars <- readRDS(file = paste0(file.name,"_bg_model_pars.rds"))
    }

  }else{

    #Select pre-compiled background
    thrs.data <- readRDS(list.files(path = precBack,pattern = "_th_dist.rds",full.names = T))
    bg.model.pars <- readRDS(list.files(path = precBack,pattern = "_bg_model_pars.rds",full.names = T))
    sfback <- read.delim(list.files(precBack,pattern = "sf",full.names = T),sep=" ",row.names = 1,header = F)
    #bimod.genes <- readRDS(list.files(path = precBack,pattern = "_bimod.rds",full.names = T))

    #@TODO
    #Ensure compatibility Ensembl ID (human backs) / gene names (mouse backs) - if not convert background to the query type (ID or name)
    # if(str_detect(rownames(data)[1],"^ENS")){
    #   #Ensembl IDs
    #   if(org == "mouse"){
    #     #If organism is mouse & gene are Ensembl ID in the data - need to convert mouse backgrounds to ensembl IDs
    #     require(org.Mm.eg.db)
    #     geneNames <- unique(c(names(thrs.data),rownames(bg.model.pars)))
    #     dfCoresp <- select(org.Mm.eg.db, keys = geneNames, keytype = 'SYMBOL',columns = c('SYMBOL', 'ENSEMBL'))
    #     #Remove duplicates & NAs
    #     dfCoresp <- dfCoresp[!is.na(dfCoresp$SYMBOL),]
    #     dfCoresp <- dfCoresp[!is.na(dfCoresp$ENSEMBL),]
    #     dfCoresp <- dfCoresp[!duplicated(dfCoresp[,c("SYMBOL")]),]
    #     dfCoresp <- dfCoresp[!duplicated(dfCoresp[,c("ENSEMBL")]),]

    #     #Update thrs.data list
    #     namesthrs <- names(thrs.data)
    #     thrs.data.updated <- sapply(seq(1,length(thrs.data)),function(x){
    #         coresp <- dfCoresp$ENSEMBL[which(dfCoresp$SYMBOL == namesthrs[x])]
    #         if(length(coresp)>0){
    #           tmp <- list(thrs.data[[x]])
    #           names(tmp) <- coresp
    #           return(tmp)
    #         }
    #       })
    #     thrs.data.updated <- thrs.data.updated[which(!sapply(thrs.data.updated, is.null))]
    #     thrs.data.updated <- sapply(seq(1,length(thrs.data.updated)),function(x){
    #       return(thrs.data.updated[[x]])
    #     })
    #     thrs.data <- thrs.data.updated
    #     rm(thrs.data.updated)

    #     #Update bg.model.pars
    #     idbg <- rownames(bg.model.pars)
    #     bg.model.pars.updated <- do.call("rbind",lapply(X = 1:nrow(bg.model.pars),FUN = function(i){
    #       coresp <- dfCoresp$ENSEMBL[which(dfCoresp$SYMBOL == idbg[i])]
    #       if(length(coresp)>0){
    #         return(c(coresp,as.numeric(bg.model.pars[i,])))
    #       }
    #     }))
    #     bg.model.pars <- bg.model.pars.updated[,2:4]
    #     bg.model.pars <- apply(bg.model.pars,2,as.numeric)
    #     rownames(bg.model.pars) <- as.character(bg.model.pars.updated[,1])
    #     colnames(bg.model.pars) <- c("theta","(Intercept)","log_umi")
    #     rm(bg.model.pars.updated)

    #     #Update sf backs
    #     namessf <- rownames(sfback)
    #     namessf.updated <- do.call("rbind",lapply(X = 1:nrow(sfback),FUN = function(i){
    #       coresp <- dfCoresp$ENSEMBL[which(dfCoresp$SYMBOL == namessf[i])]
    #       if(length(coresp)>0){
    #         return(c(coresp,sfback$V2[i]))
    #       }
    #     }))
    #     sfback <- data.frame("V2"=as.numeric(namessf.updated[,2]))
    #     rownames(sfback) <- namessf.updated[,1]
    #     rm(namessf.updated)
    #   }
    # }else{
    #   #Gene names
    #   if(org == "human"){
    #     #If organism is human & gene are names in the data - need to convert human backgrounds to gene names
    #     require(org.Hs.eg.db)
    #     ensmblIDs <- unique(c(names(thrs.data),rownames(bg.model.pars)))
    #     dfCoresp <- select(org.Hs.eg.db, keys = ensmblIDs, keytype = 'ENSEMBL',columns = c('SYMBOL', 'ENSEMBL'))
    #     #Remove duplicates & NAs
    #     dfCoresp <- dfCoresp[!is.na(dfCoresp$SYMBOL),]
    #     dfCoresp <- dfCoresp[!is.na(dfCoresp$ENSEMBL),]
    #     dfCoresp <- dfCoresp[!duplicated(dfCoresp[,c("SYMBOL")]),]
    #     dfCoresp <- dfCoresp[!duplicated(dfCoresp[,c("ENSEMBL")]),]

    #     #Update thrs.data list
    #     namesthrs <- names(thrs.data)
    #     thrs.data.updated <- sapply(seq(1,length(thrs.data)),function(x){
    #         coresp <- dfCoresp$SYMBOL[which(dfCoresp$ENSEMBL == namesthrs[x])]
    #         if(length(coresp)>0){
    #           tmp <- list(thrs.data[[x]])
    #           names(tmp) <- coresp
    #           return(tmp)
    #         }
    #       })
    #     thrs.data.updated <- thrs.data.updated[which(!sapply(thrs.data.updated, is.null))]
    #     thrs.data.updated <- sapply(seq(1,length(thrs.data.updated)),function(x){
    #       return(thrs.data.updated[[x]])
    #     })
    #     thrs.data <- thrs.data.updated
    #     rm(thrs.data.updated)

    #     #Update bg.model.pars
    #     namesbg <- rownames(bg.model.pars)
    #     bg.model.pars.updated <- do.call("rbind",lapply(X = 1:nrow(bg.model.pars),FUN = function(i){
    #       coresp <- dfCoresp$SYMBOL[which(dfCoresp$ENSEMBL == namesbg[i])]
    #       if(length(coresp)>0){
    #         return(c(coresp,as.numeric(bg.model.pars[i,])))
    #       }
    #     }))
    #     bg.model.pars <- bg.model.pars.updated[,2:4]
    #     bg.model.pars <- apply(bg.model.pars,2,as.numeric)
    #     rownames(bg.model.pars) <- as.character(bg.model.pars.updated[,1])
    #     colnames(bg.model.pars) <- c("theta","(Intercept)","log_umi")
    #     rm(bg.model.pars.updated)

    #     #Update sf backs
    #     namessf <- rownames(sfback)
    #     namessf.updated <- do.call("rbind",lapply(X = 1:nrow(sfback),FUN = function(i){
    #       coresp <- dfCoresp$SYMBOL[which(dfCoresp$ENSEMBL == namessf[i])]
    #       if(length(coresp)>0){
    #         return(c(coresp,sfback$V2[i]))
    #       }
    #     }))
    #     sfback <- data.frame("V2"=as.numeric(namessf.updated[,2]))
    #     rownames(sfback) <- namessf.updated[,1]
    #     rm(namessf.updated)
    #   }
    # }
  }

  #############################################################################

  ################### GENE QUANTIFICATION ####################

  # Keep cell type/subtype query with at least 10 cells
  query.cell.type.freq <- as.data.frame(table(metadata$cell.type))
  query.cell.type.to.keep <- as.character(query.cell.type.freq[which(query.cell.type.freq$Freq >= 10),"Var1"])
  metadataQuery <- metadata[which(metadata$cell.type %in% query.cell.type.to.keep),]
  dataQuery <- data[,as.character(metadataQuery$cell.id)]

  if(length(query.cell.type.to.keep) == 0){
    stop("Data contains no cell (sub)type with at least 10 cells")
  }else{
    #Discretize cell type queries gene expression
    if(discretize){
      #Fix seed
      set.seed(fixseed)
      cat("Discretising data ...")
      lapply(X = query.cell.type.to.keep,FUN = function(pop){
        if(!file.exists(paste0("./Query_",pop,"/",pop,"_DiscretisedData.rds"))){
          cat("Celltype : ",pop,"\n")
          pop.meta <- metadataQuery[which(metadataQuery$cell.type == pop),]
          pop.data <- dataQuery[,as.character(pop.meta$cell.id)]
          #If no pre-compiled background is used
          if(is.na(precBack)){
            #Ensure background exists
            if(file.exists(paste0(dir,"/Ref_",file.name,"/",file.name,"_th_dist.rds"))){
              pop.scaled <- suppressWarnings(PrepQueryData(query_dat = pop.data,
                                                           Ref_data = paste0(dir,"/Ref_",file.name),
                                                           dir=paste0(dir,"/Ref_",file.name,"/Query_",pop),
                                                           ncores=ncores,
                                                           file.id = pop,
                                                           saveData = T,
                                                           return.data = T,
                                                           bg.model.pars = bg.model.pars))
              pop.discrete <- Discretise_data(query_ScaledData = pop.scaled,file.id = pop,Ref_data = "./",saveData = T,return.data = F)
            }else{
              stop("No background reference found")
            }
          }else{
            #@TODO manage if background path wrong
            setwd(dir)
            pop.scaled <- suppressWarnings(PrepQueryData(query_dat = pop.data,
                                                         Ref_data = precBack,
                                                         dir=paste0(dir,"/Query_",pop),
                                                         sfback=sfback,
                                                         ncores=ncores,
                                                         file.id = pop,
                                                         saveData = T,
                                                         return.data = T,
                                                         bg.model.pars = bg.model.pars))
            pop.discrete <- Discretise_data(query_ScaledData = pop.scaled,file.id = pop,Ref_data = precBack,thsback=thrs.data,saveData = T,return.data = F)
          }
        }
      })
    }

    #Significant levels of gene expression
    if(sig.frames){
       if(!is.na(precBack)){
          setwd(dir)
       }
      #If cell type discretized
      if(file.exists(paste0("./Query_",query.cell.type.to.keep[1],"/",query.cell.type.to.keep[1],"_DiscretisedData.rds"))){
        #Gene level by population
        dir.create("sig_frames")
        #Key identity genes
        dir.create("identity_genes")
        pop.sig.frames <- lapply(X = query.cell.type.to.keep,FUN = function(pop){
          cat("Cell type : ",pop,"\n")
          discrete.data <- readRDS(file = paste0("./Query_",pop,"/",pop,"_DiscretisedData.rds"))
          discrete.mat <- discrete.data$discretised_data
          rm(discrete.data)
          invisible(gc())
          sig.frame <- Create.sig.frame(discrete.mat = discrete.mat)
          sig.frame.expr.level <- discrete.expr.levels(tbl = sig.frame[[1]], zs = sig.frame[[2]])
          #Save
          write.table(x = sig.frame.expr.level,file = paste0("./sig_frames/",pop,"_sig_tbl.txt"),sep = "\t",row.names = F,quote = F)
          #Identity genes = high level of expression & medium (if biomodal in background) sig.frame.expr.level
          #dfkig <- get.identity.genes(tbl = sig.frame.expr.level, bimodg = bimod.genes)
          #Save
          #write.table(x = dfkig,file = paste0("./identity_genes/",pop,"_identity_genes.txt"),sep = "\t",col.names = T,row.names = F,quote = F)
        })
      }
    }
  }

  ##############################################################

  return(NULL)
}
