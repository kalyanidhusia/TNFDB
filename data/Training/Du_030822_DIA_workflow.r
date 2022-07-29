

#### SWITCH TO R TERMINAL AND CREATE PROJECT DIRECTORY (LINUX COMMANDS)
## cd home/uams/analysis/
## mkdir PI_DATE

# create project directory in /home/uams/analysis/PI_DATE
# mkdir home/uams/analysis/PI_DATE

setwd("/home/uams/analysis/Du_030822_DIA")


## PASTE IN R TERMINAL
## copy project files to working directory in "files" folder
# gsutil -m cp -r "gs://uams-output/projects_2022/Du_030822_DIA" .

## LOAD FUNCTIONS
# source("/home/uams/uams-rscripts/current/functions25.R")
source("./R/functions26_CURRENT_031722-85.R")


########################################################
##  DIA ANALYSIS
########################################################
pro.file      <- "./input_files/Samples Report of Du_030822.csv"
meta.file     <- "./input_files/Du_030822_metafile_DIA.csv"
contrast.file <- NULL
pipe          <- "DIA"
enrich        <- "protein"
ilab          <- "NDu_030822_DIA"
files=list(pro.file=pro.file,meta.file=meta.file,contrast.file=contrast.file)


## CONTRAST FILES
sink(file="input_files/kidney_contrasts.txt")
cat("SLCKO_Kidney_vs_Control_Kidney = SLCKO_Kidney - Control_Kidney",sep="\n")
sink()
sink(file="input_files/brain_contrasts.txt")
cat("SLCKO_Brain_vs_Control_Brain = SLCKO_Brain - Control_Brain",sep="\n")
sink()
sink(file="input_files/intestine_contrasts.txt")
cat("SLCKO_intestine_vs_Control_intestine = SLCKO_Intestine - Control_Intestine",sep="\n")
sink()



## IMPORT SAMPLES REPORT
ext  <- extract_data(file=pro.file, sampleIDs=NULL, pipe=pipe, enrich=enrich)
head(ext$data)

## IMPORT META DATA
meta <- make_targets(file=meta.file, sampleIDs=colnames(ext$data), pipe=pipe,enrich=enrich)
meta$sample_type <- gsub(".*_","",meta$group)
head(meta)






########################################################
##    KIDNEY ANALYSIS     
########################################################

## SUBSET METADATA (KIDNEY)
sub <- subset_targets(targets=meta,factor="sample_type",rm.vals=c("Pool","Brain","Intestine"))
sub$targets

## FILTER / NORMALIZE DATA
norm <-process_data(data=ext$data, targets=sub$targets, group="group", min.reps=2, min.grps=1)

dir=file.path("protein_analysis",unique(norm$targets$sample_type),"01_quality_control");dir
pn <- make_proteinorm_report(normList=norm$normList, groups=norm$targets$group, 
                             batch=NULL, sampleLabels=NULL, legend=TRUE, enrich=enrich,
                             dir=dir, file=NULL, save=TRUE, keep.png=TRUE)

norm.meth="cycloess"
qc <- make_qc_report(normList=norm$normList, norm.meth=norm.meth,
                     groups=norm$targets$group, batch=NULL, sampleLabels=NULL, 
                     stdize=TRUE, top=500, dims=c(1,2), cex.dot=1, cex.names=1,
                     clust.metric="euclidean", clust.meth="complete",
                     xlim=NULL,ylim=NULL, enrich=enrich, dir=dir,
                     file=NULL, save=TRUE, keep.png=TRUE)


# new.dir="protein_analysis/new78"
# movelist=file.path("protein_analysis",c("01_quality_control","02_diff_expression",
#            list.files("protein_analysis",pattern="Results.xlsx")));movelist
# if(!dir.exists(new.dir)){dir.create(new.dir,recursive = T)}
# sapply(movelist, function(x){ 
#   file.copy(from=x, to=new.dir, overwrite = T, recursive = T)
#   if(all(list.files(file.path(new.dir,basename(x)),recursive=T)%in%list.files(x,recursive=TRUE))){
#     print("same removing...");unlink(x=from.dir,recursive=T)
# }})




des  <- make_design(targets=norm$targets, group="group", factors=NULL)


contrast.file="input_files/kidney_contrasts.txt"
con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts
# con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts

## create contrast.vec and the contrasts yourself 
# contrast.vec = c(paste0(group2,"_vs_",group1,"=",group2,"-",group1));contrast.vec
# contrasts <- limma::makeContrasts(contrasts=contrast.vec, levels=des$design)
# colnames(contrasts) <- gsub('=.*','',colnames(contrasts))
# con <- list(contrasts=contrasts, contrast.vec=contrast.vec)
# con$contrasts


norm.meth
dir2=file.path("protein_analysis",unique(norm$targets$sample_type),"02_diff_expression");dir2
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
                          dir        = dir2,
                          save       = TRUE,
                          ilab       = ilab  ## PI_DATE
)


wb<-openxlsx::createWorkbook()
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.meth, 
                            min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                            pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)
filename<-paste0(lim$param$Value$ilab,"_Kidney_Results.xlsx");filename
wb.file <- file.path(dirname(lim$param$Value$dir),filename); wb.file
openxlsx::saveWorkbook(wb=wb,file=wb.file, overwrite=TRUE,returnValue=TRUE)
openxlsx::saveWorkbook(wb=wb,file=file.path("logs",filename), overwrite=TRUE,returnValue=TRUE)


## SAVE RDA
if(!dir.exists("rdata")){dir.create("rdata")}
save(ext, meta, sub, norm, norm.meth, pn, qc, des, con, lim, file=file.path("rdata","lim.kidney.rda"))

cbind(lim$sum.dtp,lim$sum.dt)
#        SLCKO_Kidney_vs_Control_Kidney SLCKO_Kidney_vs_Control_Kidney
# Down                               39                              1
# NotSig                           4526                           4601
# Up                                 37                              0





########################################################
##    BRAIN ANALYSIS     
########################################################

## SUBSET METADATA (BRAIN)
sub <- subset_targets(targets=meta,factor="sample_type",rm.vals=c("Pool","Kidney","Intestine"))
sub$targets

## FILTER / NORMALIZE DATA
norm <-process_data(data=ext$data, targets=sub$targets, group="group", min.reps=2, min.grps=1)

dir=file.path("protein_analysis",unique(norm$targets$sample_type),"01_quality_control");dir
pn <- make_proteinorm_report(normList=norm$normList, groups=norm$targets$group, 
                             batch=NULL, sampleLabels=NULL, legend=TRUE, enrich=enrich,
                             dir=dir, file=NULL, save=TRUE, keep.png=TRUE)

norm.meth="cycloess"
qc <- make_qc_report(normList=norm$normList, norm.meth=norm.meth,
                     groups=norm$targets$group, batch=NULL, sampleLabels=NULL, 
                     stdize=TRUE, top=500, dims=c(1,2), cex.dot=1, cex.names=1,
                     clust.metric="euclidean", clust.meth="complete",
                     xlim=NULL,ylim=NULL, enrich=enrich, dir=dir,
                     file=NULL, save=TRUE, keep.png=TRUE)


# new.dir="protein_analysis/new78"
# movelist=file.path("protein_analysis",c("01_quality_control","02_diff_expression",
#            list.files("protein_analysis",pattern="Results.xlsx")));movelist
# if(!dir.exists(new.dir)){dir.create(new.dir,recursive = T)}
# sapply(movelist, function(x){ 
#   file.copy(from=x, to=new.dir, overwrite = T, recursive = T)
#   if(all(list.files(file.path(new.dir,basename(x)),recursive=T)%in%list.files(x,recursive=TRUE))){
#     print("same removing...");unlink(x=from.dir,recursive=T)
# }})




des  <- make_design(targets=norm$targets, group="group", factors=NULL)


contrast.file="input_files/brain_contrasts.txt"
con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts
# con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts

## create contrast.vec and the contrasts yourself 
# contrast.vec = c(paste0(group2,"_vs_",group1,"=",group2,"-",group1));contrast.vec
# contrasts <- limma::makeContrasts(contrasts=contrast.vec, levels=des$design)
# colnames(contrasts) <- gsub('=.*','',colnames(contrasts))
# con <- list(contrasts=contrasts, contrast.vec=contrast.vec)
# con$contrasts


norm.meth
dir2=file.path("protein_analysis",unique(norm$targets$sample_type),"02_diff_expression");dir2
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
                          dir        = dir2,
                          save       = TRUE,
                          ilab       = ilab  ## PI_DATE
)


wb<-openxlsx::createWorkbook()
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.meth, 
                  min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                  pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)
filename<-paste0(lim$param$Value$ilab,"_Brain_Results.xlsx");filename
wb.file <- file.path(dirname(lim$param$Value$dir),filename); wb.file
openxlsx::saveWorkbook(wb=wb,file=wb.file, overwrite=TRUE,returnValue=TRUE)
openxlsx::saveWorkbook(wb=wb,file=file.path("logs",filename), overwrite=TRUE,returnValue=TRUE)


## SAVE RDA
if(!dir.exists("rdata")){dir.create("rdata")}
save(ext, meta, sub, norm, norm.meth, pn, qc, des, con, lim, file=file.path("rdata","lim.brain.rda"))

cbind(lim$sum.dtp,lim$sum.dt)
#        SLCKO_Brain_vs_Control_Brain SLCKO_Brain_vs_Control_Brain
# Down                             35                            0
# NotSig                         4538                         4593
# Up                               20                            0






########################################################
##    INTESTINE ANALYSIS     
########################################################

## SUBSET METADATA (INTESTINE)
sub <- subset_targets(targets=meta,factor="sample_type",rm.vals=c("Pool","Kidney","Brain"))
sub$targets

## FILTER / NORMALIZE DATA
norm <-process_data(data=ext$data, targets=sub$targets, group="group", min.reps=2, min.grps=1)

dir=file.path("protein_analysis",unique(norm$targets$sample_type),"01_quality_control");dir
pn <- make_proteinorm_report(normList=norm$normList, groups=norm$targets$group, 
                             batch=NULL, sampleLabels=NULL, legend=TRUE, enrich=enrich,
                             dir=dir, file=NULL, save=TRUE, keep.png=TRUE)

norm.meth="cycloess"
qc <- make_qc_report(normList=norm$normList, norm.meth=norm.meth,
                     groups=norm$targets$group, batch=NULL, sampleLabels=NULL, 
                     stdize=TRUE, top=500, dims=c(1,2), cex.dot=1, cex.names=1,
                     clust.metric="euclidean", clust.meth="complete",
                     xlim=NULL,ylim=NULL, enrich=enrich, dir=dir,
                     file=NULL, save=TRUE, keep.png=TRUE)






des  <- make_design(targets=norm$targets, group="group", factors=NULL)


contrast.file="input_files/intestine_contrasts.txt"
con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts
# con  <- make_contrasts(file=contrast.file, design=des$design);con$contrasts

## create contrast.vec and the contrasts yourself 
# contrast.vec = c(paste0(group2,"_vs_",group1,"=",group2,"-",group1));contrast.vec
# contrasts <- limma::makeContrasts(contrasts=contrast.vec, levels=des$design)
# colnames(contrasts) <- gsub('=.*','',colnames(contrasts))
# con <- list(contrasts=contrasts, contrast.vec=contrast.vec)
# con$contrasts


norm.meth
dir2=file.path("protein_analysis",unique(norm$targets$sample_type),"02_diff_expression");dir2
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
                          dir        = dir2,
                          save       = TRUE,
                          ilab       = ilab  ## PI_DATE
)


wb<-openxlsx::createWorkbook()
add_limma_results(wb=wb, statList=lim$statList, annot=lim$annot,data=lim$data, norm.method=norm.meth, 
                  min.pval=lim$param$Value$min.pval, min.lfc=lim$param$Value$min.lfc, 
                  pipe=lim$param$Value$pipe, enrich=lim$param$Value$enrich)
filename<-paste0(lim$param$Value$ilab,"_Intestine_Results.xlsx");filename
wb.file <- file.path(dirname(lim$param$Value$dir),filename); wb.file
openxlsx::saveWorkbook(wb=wb,file=wb.file, overwrite=TRUE,returnValue=TRUE)
openxlsx::saveWorkbook(wb=wb,file=file.path("logs",filename), overwrite=TRUE,returnValue=TRUE)


## SAVE RDA
if(!dir.exists("rdata")){dir.create("rdata")}
save(ext, meta, sub, norm, norm.meth, pn, qc, des, con, lim, file=file.path("rdata","lim.intestine.rda"))

cbind(lim$sum.dtp,lim$sum.dt)
#        SLCKO_intestine_vs_Control_intestine SLCKO_intestine_vs_Control_intestine
# Down                                     80                                    0
# NotSig                                 4401                                 4603
# Up                                      122                                    0

























# 
# 
# 
# 
# ########## run in terminal and not console!!!!!!!!!!
# 
# # copy Scaffold file from the prot_fs folder or MaxQuant /txt folder to the protein_analysis folder 
# # sitting in /home/uams/analysis/PI_DATE
# # cp - r "/home/uams/analysis/Odle_021522_iBAQ/docs/Odle_021522_iBAQ/MaxQuant" "/home/uams/analysis/Odle_021522_iBAQ/protein_analysis/"
# # cp "/home/uams/prot_fs/projects/Balachandran_080521/Balachandran_080521.sdia" "/home/uams/analysis/Balachandran_080521/protein_analysis"
# 
# 
# # must be sitting in the directory with the "protein_analysis" folder visible 
# # Ctrl+Alt+enter to send to terminal
python3 /home/uams/zip/createZip.py Nagalo 020722
# 
# # copy .zip file to uams-enduser bucket
gsutil -m cp Nagalo_020722.zip gs://uams-enduser
# 
# # copy everything to uams-output
# # from working directory (copy the project to the Analysis folder)
# gsutil -m cp -r /home/uams/analysis/Odle_021522_iBAQ/ gs://uams-output/projects_2022/Odle_021522_iBAQ/Analysis
# 
gsutil -m cp -r /home/uams/analysis/Nagalo_020722_FFPE_DDA/ gs://uams-output/projects_2022/Nagalo_020722_FFPE_DDA/Analysis


