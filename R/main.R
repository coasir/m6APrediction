#' @title Encode DNA Sequences into Factor Data Frame
#'
#' @description This function takes a vector of DNA strings (e.g., 5-mers) and converts them
#' into a data frame where each nucleotide position is represented as a factor.
#' This format is necessary for many machine learning models in R.
#'
#' @param dna_strings A character vector where each element is a DNA sequence
#'   (all sequences must be of the same length).
#'
#' @return A data frame where the number of columns equals the length of the
#'   DNA strings. Each column is named 'nt_posX' (where X is the position)
#'   and is a factor with levels 'A', 'T', 'C', 'G'.
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' @title Predict m6A Status for Multiple Samples
#'
#' @description Uses a trained machine learning model to predict the m6A modification
#' probability and final status (Positive/Negative) for a batch of samples.
#'
#' @importFrom randomForest randomForest
#' @importFrom stats predict
#' @export
#'
#' @param ml_fit The trained machine learning model object (e.g., from the 'randomForest' package).
#' @param feature_df A data frame containing all required feature columns:
#'   "gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction",
#'   "evolutionary_conservation", and "DNA_5mer".
#' @param positive_threshold The probability threshold (between 0 and 1) above which
#'   a site is classified as "Positive" (m6A). Default is 0.5.
#'
#' @return The input data frame (`feature_df`) with two new columns appended:
#'   `predicted_m6A_prob` (numeric probability) and `predicted_m6A_status` (character: "Positive" or "Negative").
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
  feature_df_processed <- feature_df
  feature_df_processed <- cbind(feature_df_processed, dna_encoding(feature_df_processed$DNA_5mer))
  feature_df_processed$RNA_type <- factor(feature_df_processed$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df_processed$RNA_region <- factor(feature_df_processed$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  pred_probabilities <- predict(ml_fit, newdata = feature_df_processed, type = "prob")
  positive_prob <- pred_probabilities[, "Positive"]
  predicted_status <- ifelse(positive_prob > positive_threshold, "Positive", "Negative")
  feature_df$predicted_m6A_prob <- positive_prob
  feature_df$predicted_m6A_status <- predicted_status
  return(feature_df)
}

#' @title Predict m6A Status for a Single Sample
#'
#' @description Predicts the m6A modification probability and status for a single
#' site given all its feature values.
#'
#' @export
#'
#' @param ml_fit The trained machine learning model object.
#' @param gc_content The GC content (numeric, 0 to 1).
#' @param RNA_type The type of RNA (character: "mRNA", "lincRNA", "lncRNA", or "pseudogene").
#' @param RNA_region The region of the RNA (character: "CDS", "intron", "3'UTR", or "5'UTR").
#' @param exon_length The length of the exon (numeric).
#' @param distance_to_junction The distance to the nearest splice junction (numeric).
#' @param evolutionary_conservation The conservation score (numeric, e.g., phastCons or phyloP).
#' @param DNA_5mer The 5-mer DNA sequence at the site (character, e.g., "GGACC").
#' @param positive_threshold The probability threshold for positive classification. Default is 0.5.
#'
#' @return A named character vector with two elements:
#'   \itemize{
#'     \item \code{predicted_m6A_prob}: The predicted m6A modification probability (numeric value as a string).
#'     \item \code{predicted_m6A_status}: The predicted status ("Positive" or "Negative").
#'   }
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  feature_df <- data.frame(
    gc_content = as.numeric(gc_content),
    RNA_type = as.character(RNA_type),
    RNA_region = as.character(RNA_region),
    exon_length = as.numeric(exon_length),
    distance_to_junction = as.numeric(distance_to_junction),
    evolutionary_conservation = as.numeric(evolutionary_conservation),
    DNA_5mer = as.character(DNA_5mer),
    stringsAsFactors = FALSE
  )
  pred_df <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )
  returned_vector <- c(
    "predicted_m6A_prob" = pred_df$predicted_m6A_prob[1],
    "predicted_m6A_status" = pred_df$predicted_m6A_status[1]
  )
  return(returned_vector)
}
