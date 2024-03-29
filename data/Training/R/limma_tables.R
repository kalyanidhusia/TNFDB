#' Write tables of limma results
#'
#' This wrapper function prepares and then saves a number of tables of results:
#' \enumerate{
#'   \item A .csv summarizing the number of direction of differentially
#'     expressed genes for each contrast.
#'   \item An Excel spreadsheet of the normalized expression and
#'     statistical results for all contrasts.
#'   \item A results .csv file with the same data as the Excel spreadsheet,
#'     but without all the formatting.
#'   \item A results .csv for importing into BiqQuery. This is the same as the
#'     other results .csv, but with all missing data covnerted to 0.
#'   \item A folder of results .csv file which contain the results for each
#'     contrast.
#' }
#'
#' @param model_results A model results object returned by
#'   \code{\link{extract_limma_DE_results}}.
#' @param annotation A dataframe containing annotation data for this analysis.
#'   In our pipeline, usually the "annot" slot from the object returned by
#'   \code{\link{extract_data}}.
#' @param ilab The ilab identifier for this project.
#' @param out_dir The directory in which to output tables. If not specified,
#'   will construct a directory based on the ilab name and type of analysis.
#' @param norm.method The method that was used to normalize the data for the
#'   statistical model being output.
#' @param enrich The type of analysis being run. Options are "protein" and "phospho".
#' @param pipe The analysis pipeline being run. Options are "DIA", "TMT",
#'   "phosphoTMT", and "LF".
#' @param overwrite Should results files be overwritten? Default is F.
#' @param contrasts_subdir The subdirectory within out_dir to write the per-contrast
#'   result .csv files. If not specified, will be "per_contrast_results".
#' @param summary_csv The filename of the csv file giving a summary of the number
#'   of DE genes per contrast. If not specified, will be "DE_summary.csv".
#' @param combined_file_csv The filename of the combined results csv file. If
#'   not specified, will be combined_results.csv.
#' @param BQ_csv The filename of the BiqQuery results csv file. If not specified,
#'   will construct a name based on the ilab identifier.
#' @param spreadsheet_xlsx The filename of the Excel spreadsheet containing the
#'   results. If not specified, will construct a name based on the ilab identifier.
#' @param add_filter Should per-column filters be added to the spreadsheet?
#'   Default is TRUE. These can sometimes cause unstable spreadsheet files,
#'   try setting to FALSE if you're having issues with the Excel output.
#'
#'
#' @return Invisibly returns a dataframe with the combined results.
#'
#' @export
#'
#' @examples
#' # No examples yet
write_limma_results <- function(model_results,
                                annotation,
                                ilab,
                                out_dir = NULL,
                                norm.method,
                                enrich = c("protein", "phospho"),
                                pipe = c("DIA", "TMT", "phosphoTMT", "LF"),
                                overwrite = F,
                                contrasts_subdir = NULL,
                                summary_csv = NULL,
                                combined_file_csv = NULL,
                                BQ_csv = NULL,
                                spreadsheet_xlsx = NULL,
                                add_filter = T) {

  # Check args ----------------------------------------------------
  pipe <- rlang::arg_match(pipe)
  enrich <- rlang::arg_match(enrich)


  # Assign defaults if not overridden above
  if (is.null(out_dir)) {
    out_dir <- file.path(ilab, paste0(enrich, "_analysis"))
  }
  if (is.null(contrasts_subdir)) {
    contrasts_subdir <- "per_contrast_results"
  }
  if (is.null(summary_csv)) {
    summary_csv <- "DE_summary.csv"
  }

  if (is.null(combined_file_csv)) {
    combined_file_csv <- "combined_results.csv"
  }

  if (is.null(BQ_csv)) {
    BQ_csv <- paste0(ilab, "_results_BQ.csv")
  }

  if (is.null(spreadsheet_xlsx)) {
    spreadsheet_xlsx <- paste0(ilab, "_results.xlsx")
  }

  # Validate filenames
  for (filename in c(summary_csv, combined_file_csv, BQ_csv)) {
    validate_filename(filename, allowed_exts = "csv")
  }
  validate_filename(spreadsheet_xlsx, allowed_exts = "xlsx")

  # Make sure proteins in results have matching rows in annotation
  if (!all(unique(unlist(lapply(model_results$stats_by_contrast, rownames))) %in% rownames(annotation))) {
    cli::cli_abort(c("Not all proteins/rownames in {.arg model_results} have a matching rowname in {.arg annotation}",
                     "i" = "Are these objects from the same dataset?"))
  }

  # Check norm method
  norm.method <- rlang::arg_match(
    arg = norm.method,
    values = norm.methods,
    multiple = FALSE
  )


  # Setup -------------------------------------------------------------------
  cli::cli_rule()
  if (dir.exists(out_dir)) {
    if (overwrite) {
      cli::cli_inform("Directory {.path {out_dir}} already exists. {.arg overwrite} == {.val {overwrite}}. Overwriting files in directory.")
    } else {
      cli::cli_abort(c("Directory {.path {out_dir}} already exists",
                       "!" = "and {.arg overwrite} == {.val {overwrite}}",
                       "i" = "Rename {.arg out_dir} or set {.arg overwrite} to {.val TRUE}"))
    }
  } else {
    dir.create(out_dir, recursive = T)
  }

  statlist <- model_results$stats_by_contrast
  data <- model_results$data

  # Write summary CSV -------------------------------------------------------
  summary_output_file <- file.path(out_dir, summary_csv)
  cli::cli_inform("Writing DE summary table to {.path {summary_output_file}}")

  summary <- do.call("rbind",
                     lapply(X = names(statlist),
                            FUN = summarize_contrast_DE,
                            contrast_res_list = statlist))
  summary$pval_thresh <- model_results$pval.thresh
  summary$lfc_thresh <- model_results$lfc.thresh
  summary$p_adj_method <- model_results$adj.method



  utils::write.csv(x = summary,
                   file = summary_output_file,
                   row.names = F)

  if (!file.exists(summary_output_file)) {
    cli::cli_abort(c("Failed to write summary {.path .csv} to {.path {summary_output_file}}"))
  }


  # Write per-contrast csvs -------------------------------------------------
  per_contrast_dir <- file.path(out_dir, contrasts_subdir)
  cli::cli_inform("Writing per-contrast results {.path .csv}
                  {cli::qty(length(statlist))} file{?s} to {.path {per_contrast_dir}}")


  contrast_csv_success <- write_per_contrast_csvs(annotation_df = annotation,
                                                  data = data,
                                                  results_statlist = statlist,
                                                  output_dir = per_contrast_dir)
  if (!all(contrast_csv_success)) {
    failed <- names(contrast_csv_success)[!contrast_csv_success]
    cli::cli_abort(c("Failed to write {.path .csv} results for
                     {cli::qty(sum(!contrast_csv_success))} contrast{?s}:",
                     "!" = "{.val {failed}"))
  }


  # Write combined results csv ----------------------------------------------
  combined_output_file <- file.path(out_dir, combined_file_csv)
  cli::cli_inform("Writing combined results table to {.path {combined_output_file}}")


  # Get common row order from the first element of the stat list
  row_order <- rownames(statlist[[1]])
  # Reorder rows and rename cols for each element of the results statlist
  results_for_combine <- lapply(X = names(statlist),
                                function(x) {
                                  tmp <- statlist[[x]][row_order, ,drop = F]
                                  colnames(tmp) <- paste(colnames(tmp), x, sep = "_")
                                  tmp
                                })
  combined_output_file <- file.path(out_dir, combined_file_csv)
  combined_results <- cbind(annotation[row_order, ],
                            data[row_order, ],
                            results_for_combine)
  utils::write.csv(x = combined_results,
                   file = combined_output_file,
                   row.names = F)

  if (!file.exists(combined_output_file)) {
    cli::cli_abort(c("Failed to write combined results {.path .csv} to {.path {combined_output_file}}"))
  }

  # Write combined BiqQuery CSV ---------------------------------------------
  BQ_output_file <- file.path(out_dir, BQ_csv)
  cli::cli_inform("Writing combined results table for BiqQuery to {.path {BQ_output_file}}")
  BQ_output <- combined_results
  BQ_output[is.na(BQ_output)] <- 0
  BQ_output[BQ_output == ""] <- 0

  # colnames can't begin with numbers?
  col_beings_number <- stringr::str_detect(colnames(BQ_output), "^[:digit:]")
  if (any(col_beings_number)) {
    colnames(BQ_output)[col_beings_number] <- paste0("X", colnames(BQ_output)[col_beings_number])
  }
  utils::write.csv(BQ_output, file = BQ_output_file, row.names = F)

  if (!file.exists(BQ_output_file)) {
    cli::cli_abort(c("Failed to write BigQuery results {.path .csv} to {.path {BQ_output_file}}"))
  }


  # Write excel spreadsheet -------------------------------------------------
  excel_output_file <- file.path(out_dir, spreadsheet_xlsx)
  cli::cli_inform("Writing combined results Excel spreadsheet to {.path {excel_output_file}}")



  write_limma_excel(filename = excel_output_file,
                    statlist = statlist,
                    annotation = annotation,
                    data = data,
                    norm.method = norm.method,
                    pval.thresh = model_results$pval.thresh,
                    lfc.thresh = model_results$lfc.thresh,
                    pipe = pipe,
                    enrich = enrich,
                    add_filter = add_filter)

  if (!file.exists(excel_output_file)) {
    cli::cli_abort(c("Failed to write combined results Excel spreadsheet to {.path {excel_output_file}}"))
  }


  # Finish ------------------------------------------------------------------

  cli::cli_rule()
  cli::cli_inform(c("v" = "Success"))
  # If everything works, return combined results
  invisible(combined_results)
}



#' Summarize the number of DE proteins in a contrast
#'
#' Internal function to summarize the number of DE genes/proteins for a given
#' contrast.
#'
#' @param contrast_name The name of the contrast to summarize
#' @param contrast_res_list A list of per-contrast DE results
#'
#' @return A data frame summarizing DE for the given contrast
#'
#' @examples
#' # No examples yet
summarize_contrast_DE <- function(contrast_name, contrast_res_list) {
  tmp <- contrast_res_list[[contrast_name]][,c("sig.PVal", "sig.FDR")]

  data.frame(cbind(contrast = contrast_name,
                   type = c("down", "nonsig", "up"),
                   rbind(colSums(tmp == -1, na.rm = T),
                         colSums(tmp == 0, na.rm = T),
                         colSums(tmp == 1, na.rm = T))))

}


#' Write per-contrast results csvs
#'
#' Internal function used to write per-contrast statistical results to a .csv
#' file.
#'
#' @param annotation_df A dataframe of annotation data for each gene/protein.
#' @param data A dataframe of average expression data for each sample.
#' @param results_statlist A list of per-contrast DE results.
#' @param output_dir The directory in which to save the per-contrast csv files.
#'
#' @return A logical vector indicting whether each contrast file was
#'   successfully written.
#'
#' @examples
#' # No examples yet
write_per_contrast_csvs <- function(annotation_df,
                                    data,
                                    results_statlist,
                                    output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = T)
  }

  # make list combining
  # annotation, intensity data, and statistical results
  # reordering everything by rownames to match and subset
  # Note: technically, these could be in different orders across contrasts,
  # If they differed across contrasts.
  per_contrast_results <- lapply(X = results_statlist,
                                 FUN = function(x) cbind(annotation_df[rownames(x), ],
                                                         data[rownames(x), ],
                                                         x))
  # Make corresponding filenames
  filenames <- file.path(output_dir, paste0(names(per_contrast_results), ".csv"))
  # write
  tmp <- lapply(seq_along(per_contrast_results),
                function(x)
                  {
                  utils::write.csv(per_contrast_results[[x]],
                                   file = filenames[x],
                                   row.names = F)
                  }
                )

  # check
  write_success <- file.exists(filenames)
  names(write_success) <- names(per_contrast_results)

  write_success
}


#' Write an .xlsx file of limma results
#'
#' Internal function for creating a nicely formatted Excel spreadsheet of
#' differential expression results.
#'
#' @param filename The filename of the Excel spreadsheet to be saved.
#' @param statlist A list of the per-contrast statistical results.
#' @param annotation A dataframe of annotation date for each protein.
#' @param data A dataframe containing the average expression data for each sample.
#' @param norm.method The method that was used to normalize the data for the
#'   statistical model being output.
#' @param pval.thresh The p-value threshold that was used to determine significance.
#' @param lfc.thresh The logFC threshold that was used to determine significance.
#' @param pipe The analysis pipeline being run. Options are "DIA", "TMT",
#'   "phosphoTMT", and "LF".
#' @param enrich The type of analysis being run. Options are "protein" and "phospho".
#' @param add_filter Should per-column filters be added to the spreadsheet?
#'
#' @return Invisibly returns a list, where the first element is the filename
#'   of the saved Excel spreadsheet and the second element is the openxlsx
#'   workbook object.
#'
#' @examples
#' # No examples yet
write_limma_excel <- function(filename, statlist, annotation, data, norm.method,
                              pval.thresh, lfc.thresh, pipe, enrich, add_filter) {


  # Maybe some argument processing
  normName <- names(norm.methods)[grep(norm.method,norm.methods)]


  # initialize workbook
  wb<-openxlsx::createWorkbook()

  # Set up annotation columns -----------------------------------------------
  # Can I streamline this? And how much does this need to be here?
  # Are these the same column names we enforced when making the annotation
  # in the extract data function? # excvlusivity is removed
  if(pipe == "DIA" & enrich == "protein") {
    annotCols <- c("id", "Protein.Name","Accession.Number","Molecular.Weight", "Protein.Group.Score",
                   "Identified.Peptide.Count",
                   "UniprotID", "Gene_name","Description")
    newNames <- stringr::str_replace_all(annotCols, "[._]", " ") %>%
      stringr::str_replace("UniprotID", "UniProt ID")
    annot.title <- "Protein Annotation"
    data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Intensities")
    sheetName <- "Protein Results"

  } else if (pipe == "TMT" & enrich == "protein") {

    annotCols <- c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",
                   "Gene_name","Description")
    newNames <- stringr::str_replace_all(annotCols, "[._]", " ") %>%
      stringr::str_replace("UniprotID", "UniProt ID")
    annot.title <- "Protein Annotation"
    data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Exclusive MS1 Intensities")
    sheetName <- "Protein Results"

  } else if (pipe == "LF" & enrich == "protein") {

    annotCols <- c("id", "Fasta.headers","Majority.protein.IDs", "Score", "UniprotID",
                   "Gene_name","Description")
    newNames <- stringr::str_replace_all(annotCols, "[._]", " ") %>%
      stringr::str_replace("UniprotID", "UniProt ID")
    annot.title <- "Protein Annotation"
    data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Intensities")
    sheetName <- "Protein Results"

  } else if (pipe == "phosphoTMT" & enrich == "protein") {

    annotCols <- c("id", "Fasta.headers", "Majority.protein.IDs", "Phospho..STY..site.IDs", "Score",
                   "UniprotID","Gene_name","Description") %>%
      stringr::str_replace("UniprotID", "Uniprot ID")
    newNames <- stringr::str_replace(annotCols, "\\.\\.S", " (S") %>%
      stringr::str_replace("Y\\.\\.", "Y) ") %>%
      stringr::str_replace_all("[._]", " ") %>%
      stringr::str_replace("UniprotID", "UniProt ID")

    annot.title <- "Protein Annotation"
    data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Exclusive MS1 Intensities")
    sheetName <- "Protein Results"

  } else if (pipe == "phosphoTMT" & enrich == "phospho") {

    annotCols <- c("id","Fasta.headers", "proGroupID","UniprotID","Gene_name",
                   "Description", "PEP","Score", "Localization.prob", "Phospho..STY..Probabilities",
                   "Flanking","phosAAPosition","Class")
    newNames <- stringr::str_replace(annotCols, "\\.\\.S", " (S") %>%
      stringr::str_replace("Y\\.\\.", "Y) ") %>%
      stringr::str_replace_all("[._]", " ") %>%
      stringr::str_replace("UniprotID", "UniProt ID")

    annot.title <- "Protein/Phospho Annotation"
    data.title <- paste("Log2", ifelse(normName == "Log2", "", normName), "Normalized Exclusive MS1 Intensities")
    sheetName <- "Phospho Results"
  } else {
    cli::cli_abort("Invalid combination of {.arg pipe} and {.arg enrich}")
  }


  # Initialize worksheet ----------------------------------------------------
  openxlsx::addWorksheet(wb, sheetName = sheetName, tabColour = "#FFE100")
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 3, heights = 30) # increase row height for the title row
  openxlsx::setRowHeights(wb, sheet = sheetName, rows = 4, heights = 20) # increase row height for the header row
  title_row <- 3

  # Add annotation columns -----------------------------------------------------

  # Subset and rename
  annot <- annotation[, annotCols]
  colnames(annot) <- newNames
  annot <- make_excel_hyperlinks(data = annot,
                                 url.col = "UniProt ID",
                                 url = "https://www.uniprot.org/uniprot/")

  ## merge annotation columns
  annot_col_start <- 1
  annot_col_end <- ncol(annot)

  # Make title row
  openxlsx::mergeCells(wb, sheet = sheetName,
                       cols = annot_col_start:annot_col_end,
                       rows = title_row)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = annot.title,
                      startCol = annot_col_start, startRow = title_row,
                      colNames = FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill = "#636262",
                                      halign = "center", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "white", fontSize = 14,
                                      numFmt = "TEXT")
  openxlsx::addStyle(wb = wb, sheet = sheetName,
                     style = mergeStyle, rows = 3, cols = 1, stack = TRUE)

  ## write annotation data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#d8d8d8",
                                    halign = "left", valign = "center",
                                    textDecoration = "bold",
                                    fontColour = "black", fontSize=11,
                                    numFmt="TEXT")
  openxlsx::addStyle(wb, sheet = sheetName, style = colStyle,
                     rows = title_row + 1,
                     cols = annot_col_start:annot_col_end, stack = TRUE)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = annot,
                      startCol = annot_col_start,
                      startRow = title_row + 1,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders = "columns",
                      keepNA = TRUE,
                      na.string = "NA",
                      sep = "\t")

  ## add hyperlink style (blue font, underlined)
  hyperlinkStyle <- openxlsx::createStyle(fontColour = "#0000FF",
                                          halign = "left",
                                          valign = "center",
                                          textDecoration = "underline")
  cols <- grep("UniProt ID", colnames(annot))
  openxlsx::addStyle(wb, sheet = sheetName,
                     style = hyperlinkStyle, cols = cols,
                     rows = (title_row + 2):(nrow(annot) + 4),
                     gridExpand = TRUE, stack = TRUE)




  # Add data columns --------------------------------------------------------
  ## start/stop column positions for data
  data_col_start <- annot_col_end + 1
  data_col_end <- annot_col_end + ncol(data)

  ## Set up title row
  openxlsx::mergeCells(wb, sheet = sheetName,
                       cols = data_col_start:data_col_end,
                       rows = title_row)
  openxlsx::writeData(wb, sheet = sheetName,
                      x = data.title,
                      startCol = data_col_start,
                      startRow = title_row, colNames = FALSE)
  mergeStyle <- openxlsx::createStyle(fgFill = "#516285",
                                      halign = "center", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "white",
                                      fontSize = 14, numFmt = "TEXT")
  openxlsx::addStyle(wb, sheet = sheetName,
                     style = mergeStyle,
                     rows = title_row, cols = data_col_start, stack = TRUE)

  ## write data to worksheet
  colStyle <- openxlsx::createStyle(fgFill = "#cdd4e1",
                                    halign = "left", valign = "center",
                                    textDecoration = "bold",
                                    fontColour = "black", fontSize = 11,
                                    numFmt = "TEXT")
  openxlsx::addStyle(wb, sheet = sheetName, style = colStyle,
                     rows = title_row + 1, cols = data_col_start:data_col_end, stack=TRUE)
  openxlsx::writeData(wb, sheet = sheetName, x = data,
                      startCol = data_col_start,
                      startRow = title_row + 1,
                      colNames = TRUE,
                      rowNames = FALSE,
                      borders = "columns",
                      keepNA = TRUE, na.string = "NA",
                      sep="\t")


  # Add stats columns -------------------------------------------------------
  # colors for stat columns
  if (length(names(statlist)) > length(binfcolors)) {
    numColors<- ceiling(length(names(statlist))/length(binfcolors)) + 1
    colors2 <- rep(binfcolors, numColors)
    lightcolors2 <- rep(lightbinfcolors, numColors)
  } else {
    if (length(names(statlist)) <= length(binfcolors)) {
      colors2 <- binfcolors
      lightcolors2 <- lightbinfcolors
    }
  }
  stopifnot(length(colors2) >= length(names(statlist)))


  ## write stat data to worksheet
  stat.start<- NULL
  stat.end  <- NULL

  for(i in base::seq_along(names(statlist))) {

    stats <- statlist[[i]]
    comparison <- names(statlist)[i]
    title_color <- colors2[i]
    header_color <- lightcolors2[i]

    ## stats start
    if(is.null(stat.start)){
      stat.start = data_col_end + 1
    } else {
      if(!is.null(stat.start)){
        stat.start = stat.end + 1}
    }

    ## stats end
    if (is.null(stat.end)) {
      stat.end = data_col_end + ncol(stats)
    } else {
      if (!is.null(stat.end)) {
        stat.end = stat.end + ncol(stats)
      }
    }

    # set up title row
    openxlsx::mergeCells(wb, sheet = sheetName,
                         cols = stat.start:stat.end, rows = title_row)
    openxlsx::writeData(wb, sheet = sheetName,
                        x = comparison,
                        startCol = stat.start, startRow = title_row, colNames = FALSE)
    mergeStyle <- openxlsx::createStyle(fgFill = title_color, halign = "center",
                                        valign = "center", textDecoration = "bold",
                                        fontColour = "white", fontSize = 14, numFmt = "TEXT")
    openxlsx::addStyle(wb, sheet = sheetName, style = mergeStyle, rows = title_row, cols = stat.start, stack=TRUE)


    ## write stats to worksheet
    colStyle <- openxlsx::createStyle(fgFill = header_color,
                                      halign = "left", valign = "center",
                                      textDecoration = "bold",
                                      fontColour = "black", fontSize = 11, numFmt = "TEXT")
    openxlsx::addStyle(wb, sheet = sheetName, style = colStyle, rows = title_row + 1,
                       cols = stat.start:stat.end, stack = TRUE)
    openxlsx::writeData(wb, sheet = sheetName, x = stats,
                        startCol = stat.start,
                        startRow = title_row + 1,
                        colNames = TRUE,
                        rowNames = FALSE,
                        borders = "columns",
                        keepNA = TRUE, na.string = "NA",
                        sep = "\t")

    # Add conditional formatting to stats section
    posStyle  <- openxlsx::createStyle(fontColour="#990000", bgFill="#FFC7CE")
    negStyle  <- openxlsx::createStyle(fontColour="#006100", bgFill="#C6EFCE")
    normStyle <- openxlsx::createStyle(fontColour="#000000", bgFill="#FFFFFF")

    fc.col    <- grep("logFC", colnames(stats)) + stat.start - 1
    fc.rule1  <- paste0(">=", lfc.thresh)
    fc.rule2  <- paste0("<=", -lfc.thresh)
    fdr.col   <- grep("adj.P.Val", colnames(stats)) + stat.start - 1
    fdr.rule  <- paste0("<=", pval.thresh)
    pval.col  <- grep("P.Value", colnames(stats)) + stat.start - 1
    pval.rule <- paste0("<=", pval.thresh)


    openxlsx::conditionalFormatting(wb, sheet = sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) + 4),
                                    type = "expression", rule = fc.rule1, style = posStyle)
    openxlsx::conditionalFormatting(wb = wb, sheet = sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) + 4),
                                    type="expression", rule = fc.rule2, style = negStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName,
                                    cols = fc.col, rows = (title_row + 2):(nrow(stats) +4),
                                    type="contains", rule="NA", style = normStyle)

    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=fdr.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=fdr.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)

    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="expression", rule=pval.rule, style=posStyle)
    openxlsx::conditionalFormatting(wb=wb, sheet=sheetName, cols=pval.col, rows=5:(nrow(stats)+4),
                                    type="contains", rule="NA", style=normStyle)


  }

  ## add global styles
  headerStyle <- openxlsx::createStyle(halign = "center", valign="center", #fontColour = "#000000",
                                       border = "TopBottomLeftRight", borderColour = c("black","black","black","black"),
                                       textDecoration = "bold")
  openxlsx::addStyle(wb, sheet = sheetName, style = headerStyle,
                     cols = 1:stat.end, rows = title_row,
                     gridExpand = TRUE, stack = TRUE)
  openxlsx::addStyle(wb, sheet = sheetName, style = headerStyle,
                     cols = 1:stat.end, rows = title_row + 1,
                     gridExpand = TRUE, stack = TRUE)
  if (add_filter) {
    openxlsx::addFilter(wb, sheet = sheetName, cols = 1:stat.end, rows = title_row + 1)
  }


  # Save workbook
  openxlsx::saveWorkbook(wb = wb,
                         file = filename,
                         overwrite = TRUE)

  invisible(list(file = filename,
                 wb = wb))
}


#' Make hyperlinks for an excel column
#'
#' Internal functions to add hyperlinks to a column in a data frame to
#' be exported to an excel file. Implicitly, works only on one column.
#'
#' @param data The data frame in which to add hyperlinks
#' @param url.col The name of the column to add hyperlinks to
#' @param url The url that will be prepended to the info in the column
#'
#' @return The original dataframe, now with hyperlinks in the desired column
#'
#' @examples
#' # No examples yet
make_excel_hyperlinks <- function(data, url.col, url) {

  ids <- data[,url.col]
  tmp <- is.na(ids)
  url2 <- paste0(url, ids)
  ids2 <- paste0("HYPERLINK(\"", url2, "\", \"", ids, "\")")
  ids2[tmp] <- NA
  data[,url.col] <- ids2
  class(data[,url.col]) <- "formula"

  return(data)
}

