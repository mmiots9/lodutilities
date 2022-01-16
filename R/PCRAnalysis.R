#' @name PCRAnalysis
#' @title PCR Analysis Tool
#' @author Matteo Miotto
#' @description This functions takes as input the results file from a PCR analysis and gives back an analysis sheet
#' @param inputfile Character string of input file, can be left NA (default) and a prompt will ask for file
#' @param outputfile Character string of output file, if left NA (default) sheet is added to input file
#' @param max_rep_diff Maximum range of replicates values
#' @param housekeeping Character string of the housekeeping gene (GADPH as default) against which all deltaCT are calculated
#' @param sheetname Character string of the output sheetname (Analysis as default)
#' @usage PCRAnalysis(inputfile = NA, outputfile = NA, max_rep_diff = 0.8, housekeeping = "GADPH", sheetname = "Analysis")
#' @returns An xlsx sheet with analysis performed (deltaCT and foldchange over housekeeping gene for each sample)

#' @importFrom xlsx read.xlsx write.xlsx
#' @importFrom dplyr select rename mutate
#' @importFrom tidyr spread
#' @importFrom magrittr %>%

#' @export
PCRAnalysis <- function(inputfile = NA, outputfile = NA, max_rep_diff = 0.8, housekeeping = "GADPH",
                        sheetname = "Analysis"){

  # Input QC
    stopifnot((typeof(max_rep_diff) == "double") & (typeof(housekeeping) == "character"))

    if (is.na(inputfile)){
      cat("Choose input file:")
      Sys.sleep(1)
      inputfile <- file.choose()
    }

    stopifnot(typeof(inputfile) == "character")

  # load file
    data <- read.xlsx(inputfile, sheetName = "Results")

  # delete first rows
    colrow <- apply(X = data, MARGIN = 1, function(x) any(is.na(x)))
    colrow <- which(!colrow)

    colnames <- data[colrow,]
    data <- data[-c(1:colrow),]
    colnames(data) <- colnames

  # Subset column of interest
    data <- data %>%
      select(`Sample Name`, `Target Name`, CT) %>%
      na.omit() %>%
      rename("Sample"=`Sample Name`, "Target"=`Target Name`)

  # stop if housekeeping not in Target
    stopifnot(housekeeping %in% data$Target)

  # change typeof CT, undetermined -> NA
    data$CT[which(data$CT == "Undetermined")] <- NA
    data$CT <- as.numeric(data$CT)

  # split on sample
    splitted <- split(data, data$Sample)


  # loop
    for (s in seq_along(splitted)) {
      df <- splitted[[s]]

      df <- df %>%
        mutate(CTn = rep(c(1:3), length(CT)/3)) %>%
        select(Target, CT, CTn) %>%
        spread(CTn, CT)

      rownames(df) <- df$Target
      df$Target <- NULL

      # substitute NA with 50
      df[is.na(df)] <- 50

      # reorder replicates
      df <- t(apply(df, 1, sort))

      df[which(df == 50)] <- NA

      # change colnames
      colnames(df) <- c(paste0("Rep", seq(1:(length(colnames(df))-1))))

      # order rows to have housekeeping first
      df <- rbind(df[which(rownames(df) == housekeeping),], df[-which(rownames(df) == housekeeping),])
      rownames(df)[1] <-  housekeeping
      df <- as.data.frame(df)

      # check for differences in replicates
      delmin <- df[,2] - df[,1] > (max_rep_diff/2)
      delmin[is.na(delmin)] <- F

      delmax <- df[,dim(df)[2]] - df[,(dim(df)[2]-1)] > (max_rep_diff/2)
      delmax[is.na(delmax)] <- F

      # subsitute with NA values that exceed replicates threshold
      df[delmin, 1] <- NA
      df[delmax, dim(df)[2]] <- NA

      # calculate values
        df <- df %>%
          mutate(mean = apply(df, 1, mean, na.rm = T), stdv = apply(df, 1, sd, na.rm = T))
        df <- df %>%
          mutate(deltaCT = df$mean - df$mean[1], foldchange = format(2^(-deltaCT), scientific=F))

      # create dataframe final

        df <- df %>%
          mutate(target = rownames(df)) %>%
          select(target, colnames(df))
        df$foldchange <- as.numeric(df$foldchange)

        rownames(df) <- NULL
        splitted[[s]] <- df

        final <- df %>%
          mutate(sample = c(names(splitted[s]), rep("", (length(df$target)-1)))) %>%
          select(sample, colnames(df))

        if (s == 1){
         newf <- final
        } else {
          newf <- cbind(newf, cbind(" " = rep(" ", length(final$sample)), final))
        }

    }

  if (is.na(outputfile)) {
    write.xlsx(newf, inputfile, sheetName = sheetname, append = T, row.names = F, showNA = F)
  } else {
    write.xlsx(newf, outputfile, sheetName = sheetname, row.names = F, showNA = F)
  }


}
