#' Encode DNA strings into a data.frame of nucleotide factors
#'
#' This function takes a character vector of equal-length DNA strings (composed of the characters “A”, “T”, “C”, and “G”) and returns a data.frame in which each column corresponds to one nucleotide position (e.g., nt_pos1, nt_pos2, …). Each cell is a factor with levels `c("A", "T", "C", "G")`.
#'
#' @param dna_strings A character vector of DNA sequences. All sequences must have the same length.
#' @return A data.frame with one row per element of `dna_strings`, and one column per nucleotide position. Each column is a factor with levels “A”, “T”, “C”, “G”.
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])
  seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict multiple new observations using a trained model
#'
#' This function uses a trained model object (`ml_fit`) to generate predictions on a new data frame (`feature_df`) of multiple observations. It enforces that certain required columns are present, recodes categorical variables, expands a 5‑mer DNA string into individual nucleotide‑position factors, and then uses `predict(..., type = "prob")` to compute the probability of the “Positive” class. It then assigns a status (“Positive” or “Negative”) based on the `positive_threshold`. The output is the original `feature_df` with two additional columns: `predicted_m6A_prob` and `predicted_m6A_status`.
#'
#' @param ml_fit A trained model object (e.g., from randomForest) that supports `predict(..., type = "prob")`.
#' @param feature_df A data.frame of new observations. Must include columns: `gc_content`, `RNA_type`, `RNA_region`, `exon_length`, `distance_to_junction`, `evolutionary_conservation`, `DNA_5mer`.
#' @param positive_threshold Numeric. (default = 0.5) A threshold above which the predicted probability is labelled “Positive”; otherwise “Negative”.
#' @return A data.frame with the same number of rows as `feature_df`, plus two new columns: `predicted_m6A_prob` (numeric) and `predicted_m6A_status` (character “Positive”/“Negative”).
#' @examples
#' ml_model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#' new_features <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                     package="m6APrediction"))
#'   result_df <- prediction_multiple(
#'     ml_fit = ml_model,
#'     feature_df = new_features,
#'     positive_threshold = 0.5
#'   )
#'   head(result_df)
#' @export
#' @import randomForest
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length",
                  "distance_to_junction", "evolutionary_conservation","DNA_5mer") %in% colnames(feature_df)))
  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))
  nn <- 5
  seq_m <- matrix(unlist(strsplit(as.character(feature_df$DNA_5mer), "")),
                  ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  feature_df <- cbind(feature_df, seq_df)
  prob_pred <- predict(ml_fit, newdata = feature_df, type = "prob")[, "Positive"]
  status_pred <- ifelse(prob_pred > positive_threshold, "Positive", "Negative")
  feature_df$predicted_m6A_prob <- prob_pred
  feature_df$predicted_m6A_status <- status_pred
  return(feature_df)
}

#' Predict a single new observation using a trained model
#'
#' This function wraps around `prediction_multiple()` to predict for a single new observation (rather than a full data frame). It builds a one‑row data.frame from the input feature values and passes it to `prediction_multiple()`, then extracts and returns a named vector with two elements: `predicted_m6A_prob` (numeric) and `predicted_m6A_status` (character).
#'
#' @param ml_fit A trained model object (e.g., from randomForest) that supports `predict(..., type="prob")`.
#' @param gc_content Numeric. The GC content of the sequence.
#' @param RNA_type Character. One of c("mRNA","lincRNA","lncRNA","pseudogene").
#' @param RNA_region Character. One of c("CDS","intron","3'UTR","5'UTR").
#' @param exon_length Numeric. The exon length of the RNA.
#' @param distance_to_junction Numeric. The distance of the site to the nearest splice junction.
#' @param evolutionary_conservation Numeric. A measure of evolutionary conservation of the site.
#' @param DNA_5mer Character. A 5‑mer DNA string (e.g., "ATCGG").
#' @param positive_threshold Numeric. (default = 0.5) A threshold above which the predicted probability is labelled “Positive”.
#' @return A named vector of length 2:
#'   * `predicted_m6A_prob` — numeric probability of the “Positive” class.
#'   * `predicted_m6A_status` — character, “Positive” or “Negative”.
#' @examples
#'
#' ml_model <- readRDS(system.file("extdata", "rf_fit.rds", package="m6APrediction"))
#'   single_pred <- prediction_single(
#'     ml_fit = ml_model,
#'     gc_content = 0.50,
#'     RNA_type = "mRNA",
#'     RNA_region = "3'UTR",
#'     exon_length = 120,
#'     distance_to_junction = 8,
#'     evolutionary_conservation = 0.7,
#'     DNA_5mer = "GATCG",
#'     positive_threshold = 0.6
#'   )
#'   single_pred
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length,
                              distance_to_junction, evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  single_row_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )
  prediction_result <- prediction_multiple(ml_fit, single_row_df, positive_threshold)
  predicted_prob <- prediction_result$predicted_m6A_prob[1]
  predicted_status <- as.character(prediction_result$predicted_m6A_status[1])
  returned_vector <- c(predicted_m6A_prob = predicted_prob,
                       predicted_m6A_status = predicted_status)
  return(returned_vector)
}

