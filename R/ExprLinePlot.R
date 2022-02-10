#' @name ExprLinePlot
#' @title Gene expression levels lineplot
#' @author Matteo Miotto
#' @description This functions create one lineplot for each gene provided from the provided dataset
#' @param dataset Gene expression dataset with a column for gene symbols (symbol) and one column for each sample in the form "neuron_stage" can be left NA (default) and a prompt will ask for file
#' @param genelist Character string of gene of interest

#' @usage ExprLinePlot(dataset = df, genelist = c("GAPDH", "CAPN6"))
#' @returns A list containing a lineplot for each gene
#'
#' @importFrom utils read.csv
#' @import ggplot2

#' @export

ExprLinePlot <- function(dataset = NULL,
                         genelist) {

  ##------------------------------------# INPUT QC #------------------------------------##
  # Check if dataset is NA
    if (is.null(dataset)) {
      datasetn <- file.choose()
      dataset <- read.csv(datasetn, header = T, sep = ",")
    }

  # Check if genelist is a chr vector
    if(!is.character(genelist)) {
      stop("genelist must be of type character")
    }

  # Check if all genes in genelist are in dataset$symbol
    if(!all(genelist %in% dataset$symbol) & length(genelist) > 1) {
      genenot <- genelist[-which(genelist %in% dataset$symbol)]

      stop(sprintf("The following genes are not in the dataset: %s",
                   paste(genenot, collapse = ", ")))
    } else if (!all(genelist %in% dataset$symbol)) {
      stop("The gene you have inserted is not in the dataset")
    }

  ##------------------------------# VARIABLES OF INTEREST #-----------------------------##

    NeuronID <- sapply(strsplit(colnames(dataset[,-1]), "_"), "[[", 1)
    Stage    <- sapply(strsplit(colnames(dataset[,-1]), "_"), "[[", 2)
    cols <- c( "X5HT3aR" = "#1C931C",
               "Lhx6" = "#7D29D6",
               "SST" = "#CE2020",
               "VIP" = "#CE2021",
               "CPN" = "#1C931C",
               "CThPN" = "#7D29D6",
               "ScPN" = "#CE2020")

    cols <- cols[unique(NeuronID)]
    plist <- list()

  ##---------------------------------# LOOP GENELIST #----------------------------------##

    for (g in genelist) {

      # SUBSETTING
        Expression <- unlist(dataset[which(dataset$symbol == g),-1])
        data.plot <- data.frame(NeuronID = NeuronID, Stage = Stage, Expression = Expression)

      # PLOT
      temp <- ggplot(data.plot) +
        geom_line(
          aes(
            x = Stage,
            y = Expression,
            color = NeuronID,
            group = NeuronID
          ),
          size = 1) +
        scale_x_discrete(expand = c(0, 0)) +
        theme_light() +
        scale_colour_manual(values = cols) +
        labs(title = sprintf("%s", g)) +
        theme(plot.title = element_text(size = 18, face = "bold"),
              axis.line = element_line(colour = "black"),
              panel.border = element_blank(),
              panel.background = element_blank())
      plist[[g]] <- temp

    }

  ##------------------------------------# RETURN #--------------------------------------##

    return(plist)
}
