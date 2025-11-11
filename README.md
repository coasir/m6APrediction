# m6APrediction: A Machine Learning Predictor for m6A RNA Modification Sites

## 💡 Overview

The m6APrediction package provides a robust, pre-trained Random Forest model for predicting the N6-methyladenosine (m6A) modification status at specific RNA sites. This tool leverages multiple genomic and sequence-based features to offer high-confidence predictions, enabling researchers to identify potential m6A modification sites in RNA transcripts with reliable probability scores and binary classification results (Positive/Negative).

## 🚀 Installation Guide

This package is designed to be installed directly from GitHub using the remotes package.

### Step 1: Install Dependencies
First, ensure you have the necessary dependencies installed:
```r
install.packages(c("devtools", "randomForest", "remotes"))
```

### Step 2: Install m6APrediction
Install the m6APrediction package from GitHub:
```r
remotes::install_github("coasir/m6APrediction")
```

## 🛠️ Usage Examples

Once installed and loaded, you can use the exported functions to predict m6A status for single or multiple sites.

### Load Package and Example Data

The package includes an example input dataset for quick testing:
```r
# Load the package
library(m6APrediction)

# Locate the example data file bundled within the package
input_file_path <- system.file("extdata", "m6A_input_example.csv", package = "m6APrediction")

# Read the data
feature_df_test <- read.csv(input_file_path, stringsAsFactors = FALSE)

# Preview the data structure
head(feature_df_test)
```

### Multiple Site Prediction

Use `prediction_multiple()` for batch prediction of multiple sites:
```r
# Predict multiple sites (using first 5 rows for demonstration)
prediction_result <- prediction_multiple(
    ml_fit = m6APrediction::ml_fit,
    feature_df = feature_df_test[1:5, ],
    positive_threshold = 0.5
)

# View results
head(prediction_result)
```

**Expected Output:**
```
  gc_content RNA_type RNA_region exon_length distance_to_junction
1  0.6616915   lncRNA     intron    0.000000            10.668885
2  0.5223881     mRNA      3'UTR   10.790348             9.330917
  evolutionary_conservation DNA_5mer predicted_m6A_prob predicted_m6A_status
1                0.01641791    GGACC              0.596             Positive
2                0.02537313    TGACC              0.008             Negative
```

### Single Site Prediction

Use `prediction_single()` for individual site analysis:
```r
# Extract a single sample for demonstration
single_sample <- feature_df_test[1, ]

# Predict single site
single_prediction_result <- prediction_single(
    ml_fit = m6APrediction::ml_fit, 
    gc_content = single_sample$gc_content,
    RNA_type = single_sample$RNA_type,
    RNA_region = single_sample$RNA_region,
    exon_length = single_sample$exon_length,
    distance_to_junction = single_sample$distance_to_junction,
    evolutionary_conservation = single_sample$evolutionary_conservation,
    DNA_5mer = single_sample$DNA_5mer
)

# View result
print(single_prediction_result)
```

**Expected Output:**
```
predicted_m6A_prob predicted_m6A_status 
           "0.596"           "Positive" 
```

## 📈 Model Performance

The pre-trained Random Forest model demonstrates excellent performance in distinguishing m6A-positive from m6A-negative sites, as illustrated by the ROC and PRC curves below:

![ROC and PRC Curves demonstrating model performance](ROC_and_PRC.png)

The model achieves high accuracy with robust sensitivity and specificity, making it suitable for reliable m6A site prediction in various RNA types and genomic contexts.

## 🔧 Model Features

The prediction model utilizes the following key features:

- **GC Content**: Nucleotide composition around the site
- **RNA Type**: mRNA, lncRNA, lincRNA, or pseudogene
- **RNA Region**: CDS, intron, 3'UTR, or 5'UTR
- **Exon Length**: Length of the containing exon
- **Distance to Junction**: Distance to nearest splice site
- **Evolutionary Conservation**: PhastCons/phyloP scores
- **DNA 5-mer**: 5-nucleotide sequence context

## 📚 References

For more information about m6A modifications and computational prediction methods, please refer to the relevant literature on epitranscriptomic modifications and machine learning approaches in RNA biology.

## 📞 Support

For questions, bug reports, or feature requests, please visit my [GitHub Issues](https://github.com/coasir/m6APrediction/issues) page.

---

**Package Information:**
- Version: 1.0.0
- License: MIT
- Author: Shengyi Gu (Shengyi.Gu23@student.xjtlu.edu.cn)
