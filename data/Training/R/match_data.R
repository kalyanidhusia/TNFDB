



##----------------------------------
##  MATCH META, DATA, ANNOT
##----------------------------------
## columns of data ordered to match rows of targets 
## (ncol data>=nrow targets,colnames data=rownames targets)
## annot ordered to match data based on rownames
## (nrow annot >= nrow data, rownames annot = rownames data)
match_data <- function(data, targets=NULL, annot=NULL, sampleLabels=NULL){
   
   obj<-list()
   
   if(is.null(targets)){
      targets <- data.frame(sampleIDs=colnames(data), group=c(rep(1,ncol(data))),row.names=colnames(data))
   }
   obj$targets<-targets
   # print(sapply(obj$targets,class))
   
   
   if(!all(rownames(obj$targets)%in%colnames(data))){
      stop("list() object cannot be created. targets (rows) and data (columns)","
            are not in the same order.")
   }
   obj$data <- data[,rownames(obj$targets)]
   
   ## if no annot df input
   if(is.null(annot)){
      annot<-data.frame(GeneIDs=rownames(obj$data),row.names=rownames(obj$data))
   }
   if(!is.null(annot)){
      if(!all(rownames(obj$data)%in%rownames(annot))){
         stop("Error: count and annot (rows)","
            are not in the same order.")
      }
      obj$annot <- data.frame(annot[rownames(obj$data), ],row.names=rownames(obj$data))
      colnames(obj$annot)<-colnames(annot)
   }
   
   if(is.null(sampleLabels)){
      sampleLabels<-colnames(obj$data)
   }
   if(any(duplicated(sampleLabels))){
      tmpLabels <- make.names(names=sampleLabels, unique=TRUE)
      first_nondups<-grep(FALSE,1:length(tmpLabels)%in%grep("\\.[0-9]*$",tmpLabels))
      tmpLabels[first_nondups]<-paste0(tmpLabels[first_nondups],".0")
      if(any(duplicated(tmpLabels))){
         tmpLabels<-make.names(names=sampleLabels, unique=TRUE)
      }
      sampleLabels<-tmpLabels      
   }
   if(length(sampleLabels) != nrow(obj$targets)){
      stop("Error: number of sampleLabels and number of rows in targets are not the same length.")
   }
   rownames(obj$targets) <- sampleLabels
   colnames(obj$data) <- sampleLabels
   
   
   return(obj)
   
}



