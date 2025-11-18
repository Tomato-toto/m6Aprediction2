# m6APrediction
## Introduction
`m6APrediction` is an R package designed for predicting m6A modification status. This tool utilizes machine learning models, specifically random forest algorithms, to analyze gene expression data and predict the presence or absence of m6A modifications. The package aims to assist bioinformatics researchers and data scientists in quickly implementing m6A predictions.
## Installation
You can install the `m6APrediction` package from GitHub using the following commands:
```r
# Uing devtools
devtools::install_github("Tomato-toto/m6APrediction")
# Or using remotes
remotes::install_github("Tomato-toto/m6APrediction")
```

Usage Example

Here is a minimal example demonstrating the usage of two prediction functions from the m6APrediction package:

```r
library(m6APrediction)

# Example input (replace with your data frame)
example_df <- data.frame(
  gc_content = c(0.54, 0.62),
  RNA_type = c("mRNA", "lncRNA"),
  RNA_region = c("exon", "intron"),
  exon_length = c(1200, 850),
  distance_to_junction = c(45, 120),
  evolutionary_conservation = c(0.82, 0.76),
  DNA_5mer = c("ATCGG", "TTCGA")
)

# Prediction using model
single_result <- prediction_single(
ml_fit = rf_fit,
gc_content = data$gc_content[1],
RNA_type = data$RNA_type[1],
RNA_region = data$RNA_region[1],
exon_length = data$exon_length[1],
distance_to_junction = data$distance_to_junction[1],
evolutionary_conservation = data$evolutionary_conservation[1],
DNA_5mer = data$DNA_5mer[1]
)
single_result

# Batch prediction
batch_results <- prediction_multiple(rf_fit, data)
head(batch_results)
```

Performance Showcase
In Practical 4, we generated ROC and PRC curves to showcase the strong performance of the model. Below are examples of these curves: ROC Curve PRC Curve
![ROC and PRC Curve](~/BIO215/Practical6/m6APrediction/ROC,PRC.png)

Author
Xinyuan.Hu
Email:Xinyuan.Hu23@student.xjtlu.edu.cn

License
This project is licensed under the MIT License.
