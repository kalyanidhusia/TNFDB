
## DIFFERENCES BETWEEN PIPELINES
# the input file / metadata file (required column fields) differ between pipelines
# required columns in metadata files differ between DIA, TMT, phosphoTMT, and LF pipelines
# phosphoTMT(protein),TMT(protein), LF(protein) use maxquant proteinGroups.txt as input
# phosphoTMT(phospho) uses maxquant phospho(STY)sites.txt file as input
# DIA(protein) uses Scafold DIA samples report.csv file as input



########################################################
##  EXAMPLE (DIA ANALYSIS)
########################################################

## LOAD FUNCTIONS
source("./R/functions26_CURRENT_031722-85.R")

pro.file      <- "./input_files/Samples Report of Du_030822.csv" ## path to samples report file or NULL choose file (Pool_# / Sample_#)
meta.file     <- "./input_files/Du_030822_metafile_DIA.csv"      ## path to metadata file or NULL to return df of info extracted from pro.file only
contrast.file <- "./input_files/contrasts.txt"                   ## path to contrast file (if file exists) NULL choose file
pipe          <- "DIA"                                           ## MS pipeline (DIA, TMT, phosphoTMT, LF)
enrich        <- "protein"                                       ## sample enrichment (protein, phospho)
ilab          <- "NDu_030822_DIA"                                ## project name (system previously used was called ilab)



## IMPORT DIA SAMPLES REPORT
ext  <- extract_data(file=pro.file, sampleIDs=NULL, pipe=pipe, enrich=enrich)
head(ext$data)


## IMPORT META DATA
meta <- make_targets(file=meta.file, sampleIDs=colnames(ext$data), pipe=pipe,enrich=enrich)
head(meta)


## SUBSET METADATA TO INCLUDE SAMPLES THAT WILL BE USED IN ThE ANALYSIS (REMOVE OUTLIERS/SPECIFIC GROUPS)
# factor=column name in targets; rm.vals=vector; values in column name to remove; e.g. 
# removes pool and empty channel samples; 
# factor="gender",rm.vals="male" (removes all male samples)
# factor="batch", rm.vals=c(2,3) (remove samples from batches 2 and 3)
sub <- subset_targets(targets=meta,factor="group", rm.vals=c("Pool","Empty_Channel"))
sub$targets


## FILTER / NORMALIZE DATA
## intensity > 0 in at least 2 replicates in 1 or more groups (group)
norm <-process_data(data=ext$data, targets=sub$targets, group="group", min.reps=2, min.grps=1)


## PROTEINORM REPORT
pn <- make_proteinorm_report(normList=norm$normList, groups=norm$targets$group, 
                             batch=NULL, sampleLabels=NULL, legend=TRUE, enrich=enrich,
                             dir=NULL, file=NULL, save=TRUE, keep.png=TRUE)

## SELECT BEST PERFORMING NORMALIZATION METHOD 
norm.methods ## 8 methods
norm.meth="cycloess" 

## QC REPORT
qc <- make_qc_report(normList=norm$normList, norm.meth=norm.meth,
                     groups=norm$targets$group, batch=NULL, sampleLabels=NULL, 
                     stdize=TRUE, top=500, dims=c(1,2), cex.dot=1, cex.names=1,
                     clust.metric="euclidean", clust.meth="complete",
                     xlim=NULL,ylim=NULL, enrich=enrich, dir=NULL,
                     file=NULL, save=TRUE, keep.png=TRUE)


## DESIGN MATRIX 
## (NO INTERCEPT MODEL (e.g. ~0+group ; ~0+group+batch)
des  <- make_design(targets=norm$targets, group="group", factors=NULL)


## MAKE CONTRASTS 

## using contrasts defined in a file
writeLines(text=c("T1_vs_CON=T1-CON", "T2_vs_CON=T2-CON"), con="./input_files/contrasts.txt")
contrast.file <- "input_files/contrasts.txt"
con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts

## via makeContrasts() function
contrast.vec <- c("T1_vs_CON=T1-CON", "T2_vs_CON=T2-CON","T3_vs_CON=T3-CON")
contrasts <- limma::makeContrasts(contrasts=contrast.vec, levels=des$design)
colnames(contrasts) <- gsub('=.*','',colnames(contrasts))
con <- list(contrasts=contrasts, contrast.vec=contrast.vec)
con$contrasts

## all pairwise contrasts 
## groups- vector,length # samples,reflects design columns
## con <- make_all_contrasts(design=des$design, groups=des$targets$group)



## LIMMA DE ANALYSIS

norm.meth
lim <- run_limma_analysis(data       = norm$normList[[norm.meth]],
                          annot      = ext$annot,
                          targets    = des$targets,
                          design     = des$design,
                          contrasts  = con$contrasts,
                          min.pval   = 0.055,
                          min.lfc    = 1,
                          adj.method = "BH",
                          paired     = FALSE,  ## TRUE if paired samples/mixed.effects model
                          pipe       = pipe,
                          enrich     = enrich,
                          dir        = FALSE,
                          save       = TRUE,
                          ilab       = ilab  ## PI_DATE
)


## SAVE LIMMA RESULTS (EXCEL) 
wb<-openxlsx::createWorkbook()
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.meth, 
                            min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                            pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)
filename<-paste0(lim$param$Value$ilab,"_Kidney_Results.xlsx");filename
wb.file <- file.path(dirname(lim$param$Value$dir),filename); wb.file
openxlsx::saveWorkbook(wb=wb,file=wb.file, overwrite=TRUE,returnValue=TRUE)



## SAVE RDA
if(!dir.exists("rdata")){dir.create("rdata")}
save(ext, meta, sub, norm, norm.meth, pn, qc, des, con, lim, file=file.path("rdata","lim.rda"))

cbind(lim$sum.dtp,lim$sum.dt)
#        SLCKO_Kidney_vs_Control_Kidney SLCKO_Kidney_vs_Control_Kidney
# Down                               39                              1
# NotSig                           4526                           4601
# Up                                 37                              0



