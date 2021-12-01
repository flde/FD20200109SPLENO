library_load <- suppressMessages(
    list(
        # Seurat 
        library(Seurat), 
        library(Matrix)
    )
)


##################
### ImportMeta ###
##################

setGeneric("ImportMeta", valueClass=c("Seurat"), def=function(so="Seurat", ...) {standardGeneric("ImportMeta")})

setMethod("ImportMeta", signature=c(so="Seurat"), definition=function(so) {
    
    if(!dir.exists("meta")) {
        
        warning("Import meta data: meta directorie does not exist")
        
    } else {
        
        message("Importing meta data ")
        
        so@meta.data <- read.csv("meta/meta.csv", row.names = 1)
        
    }
    
    return(so)
}
         
         )

#########################
### ImportAssayMatrix ###
#########################

setGeneric("ImportAssayMatrix", valueClass=c("matrix"), def=function(assay="character", slot="character", ...) {standardGeneric("ImportAssayMatrix")})

setMethod("ImportAssayMatrix", signature=c(assay="character", slot="character"), definition=function(assay, slot) {    
    
    # Initiate Seurat object 
    mtx <- Matrix::readMM(paste0("assay/", assay, "/", slot, "/matrix.mtx"))
    mtx <- as.matrix(mtx)
    mtx <- t(mtx)
    
    colnames(mtx) <- read.csv(paste0("assay/", assay, "/", slot, "/cellid.csv"))[, 1]
    rownames(mtx) <- read.csv(paste0("assay/", assay, "/", slot, "/genes.csv"))[, 1]
    
    return(mtx)  
    
}
         )

###################
### ImportAssay ###
###################

setGeneric("ImportAssay", valueClass=c("Seurat"), def=function(so="Seurat", ...) {standardGeneric("ImportAssay")})

setMethod("ImportAssay", signature=c(so="Seurat"), definition=function(so, assays=NULL, slots=NULL) {    
    
    # Start message 
    message("Importing assay data: ")

    # Assay selection 
    if(!is.null(assays)) {
        
        assays <- assays
        
    } else {
        
        assays <- list.dirs("assay", recursive=FALSE, full.names=FALSE)
        
    }



    # Import assays and slots
    for (assay in assays) {
        
        # Initialize assay slot with counts (mendatory)
        message(paste("... Importing", assay, "counts"))
        so[[assay]] <- CreateAssayObject(ImportAssayMatrix(assay, "counts"))
        
        # Slot selection
        if(!is.null(slots)) {

            slots <- slots

        } else {

            slots <- list.dirs(paste0("assay/", assay), recursive=FALSE, full.names=FALSE)

        }
        
        for(slot in slots) {
            
                if(slot!="counts") {
                    
                    message(paste("... Importing", assay, slot))
                    so <- SetAssayData(so, assay=assay, slot=slot, ImportAssayMatrix(assay, slot))
                
                }   
        } 
    }
    
    return(so)
}
         )


########################
### ImportReductions ###
########################

setGeneric("ImportReductions", valueClass=c("Seurat"), def=function(so="Seurat", ...) {standardGeneric("ImportReductions")})

setMethod("ImportReductions", signature=c(so="Seurat"), definition=function(so) {
    
    message("Importing reductions:  ")
    for(reductions in list.dirs("reductions", recursive=FALSE, full.names=FALSE)) {
        
        message(paste("... Importing", reductions))
        
        embeddings <- read.csv(paste0("reductions/", reductions, "/reduction.csv"), row.names = 1)
        embeddings <- as.matrix(embeddings)
        colnames(embeddings) <- paste0(toupper(reductions), "_", 1:ncol(embeddings))
            
        so@reductions[[reductions]] <- CreateDimReducObject(embeddings=embeddings, key=toupper(reductions))
        
    }
    
    return(so)
    
}
         
         )

##################
### dir2seurat ###
##################

setGeneric("dir2seurat", valueClass=c("logical", "NULL", "character", "Seurat"), def=function(dir="character", ...) {standardGeneric("dir2seurat")})

setMethod("dir2seurat", signature=c(dir="character"), definition=function(dir, assays=NULL, slots=NULL) {
    
    catch_dir <- getwd()
    
    message(paste("Importing Seurat object from", dir))
    setwd(dir)
    
    # Initiate Seurat object     
    so <- CreateSeuratObject(ImportAssayMatrix("RNA", "counts"), assay="RNA")
    
    # Import meta data 
    so <- ImportMeta(so)
    
    # Import assay data 
    so <- ImportAssay(so, assays, slots)
    
    # Import reductions
    so <- ImportReductions(so)
    
    setwd(catch_dir)

    return(so)
    
}
         )

################
### BuildDir ###
################

setGeneric("BuildDir", valueClass="logical", def=function(dir="character", ...) {standardGeneric("BuildDir")})

setMethod("BuildDir", signature=c(dir="character"), definition=function(dir, overwrite=FALSE) {
    
    if(dir.exists(dir) & overwrite==FALSE) {
        
        stop("Directory already exists and overwrite is set to FALSE")
        
    } else {
        
        message(paste("Creating output directory", dir))
        
        if(dir.exists(dir)) {unlink(dir, recursive=TRUE)}
        dir.create(dir, showWarnings=FALSE)
        
    }
    
}
         )

#####################
### WriteMetaData ###
#####################

setGeneric("WriteMetaData", valueClass=c("NULL", "character"), def=function(so="Seurat", dir="character", ...) {standardGeneric("WriteMetaData")})

setMethod("WriteMetaData", signature=c(so="Seurat", dir="character"), definition=function(so, dir, file="meta.csv") {
    
    if(all(dim(so@meta.data) == 0)) {
        
        warning("Writing meta data: The meta.data slot of the Seurat object is empty")
        
    } else {
        
        # Create meta data directory
        dir.create(paste0(dir, "meta"), showWarnings=FALSE)
        
        # Write meta data 
        path <- paste0(dir, "meta/", file)
        message(paste("Writing meta data to", path))
        write.csv(so@meta.data, path)
        
    }
}
         )

######################
### WriteAssayData ###
######################

setGeneric("WriteAssayData", valueClass=c("NULL", "character"), def=function(so="Seurat", dir="character") {standardGeneric("WriteAssayData")})

setMethod("WriteAssayData", signature=c(so="Seurat", dir="character"), definition=function(so, dir) {
    
    message(paste("Writing assay data for the slots:", paste(names(so@assays), collapse = ", ")))
    
    # Create assay directory 
    dir.create(paste0(dir, "assay"), showWarnings=FALSE)
    
    for(assay in names(so@assays)) {
        
        for(slot in c("counts", "data", "scale.data")) {
            
            path <- paste0(dir, "assay/", assay, "/", slot)
            
            tryCatch(
                
                {
                    message(paste("...", assay, slot))
                    
                    # Create slot directory
                    dir.create(path, recursive=TRUE, showWarnings=FALSE)
                    
                    # Write slot data
                    writeMM(t(GetAssayData(so, slot = "counts")), paste0(path, "/matrix.mtx"))
                    write.table(rownames(GetAssayData(so, slot = "counts")), row.names=FALSE, paste0(path, "/genes.csv"))
                    write.table(colnames(GetAssayData(so, slot = "counts")), row.names=FALSE, paste0(path, "/cellid.csv"))
                    write.csv(VariableFeatures(so), row.names=FALSE, paste0(path, "/variable_features.csv"))
                    
                }, 
                
            error=function(err) {
                
                # Remove slot output dir 
                unlink(path, recursive=TRUE)
                
                message(paste("Could not write", assay, slot))
                message(head(err))
                
            }
            
            )
            
        }
        
    }
    
}
         )

##########################
### WriteReductionData ###
##########################

setGeneric("WriteReductionData", valueClass=c("NULL", "character"), def=function(so="Seurat", dir="character") {standardGeneric("WriteReductionData")})

setMethod("WriteReductionData", signature=c(so="Seurat", dir="character"), definition=function(so, dir) {

    if(length(so@reductions)==0) {
        
        warning("Writing reductions: The reductions slot of the Seurat object is empty")
        
    } else {
        
        message(paste("Writing reductions:", paste(names(so@reductions), collapse = ", ")))
        
    
        for(reductions in names(so@reductions)) {
            
            message(paste("...", reductions))
            
            dir.create(paste0(dir, "reductions/", reductions), showWarnings=FALSE, recursive = TRUE)
            write.csv(so@reductions[[reductions]][[]], paste0(dir, "reductions/", reductions, "/reduction.csv"))
            
        }
    }
}
         )

##################
### seurat2dir ###
##################

setGeneric("seurat2dir", valueClass=c("character", "NULL"), def=function(so="Seurat", dir="character", ...) {standardGeneric("seurat2dir")})

setMethod("seurat2dir", signature=c(so="Seurat", dir="character"), definition=function(so, dir, overwrite=FALSE) {

    # Build root dir 
    BuildDir(dir, overwrite)
    
    # Write meta data 
    WriteMetaData(so, dir)
    
    # Write assay data 
    WriteAssayData(so, dir)
    
    # Write reduction data
    WriteReductionData(so, dir)
    
}
         )