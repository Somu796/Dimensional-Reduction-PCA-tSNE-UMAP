
library(tidyverse)
# Load necessary library
set.seed(123)  # For reproducibility

# Parameters for the dataset
n_proteins <- 100  # Number of proteins
n_samples <- 15    # Number of samples (5 Control1, 5 Control2, 5 Treatment)

# Generate random LFQ data for Control1 and Control2
control_data <- matrix(rnorm(n_proteins * 5, mean = 10, sd = 1.5), nrow = n_proteins, ncol = 5)
treatment1_data <- matrix(rnorm(n_proteins * 5, mean = 15, sd = 1.5), nrow = n_proteins, ncol = 5)

# Introduce a trend for the Treatment group (some proteins have increased or decreased expression)
treatment2_data <- control_data + matrix(rnorm(n_proteins * 5, mean = 7.5, sd = 0.5), nrow = n_proteins, ncol = 5)

# Combine the data
lfq_data <- t(cbind(control_data, treatment1_data, treatment2_data))

# Assign row and column names
colnames(lfq_data) <- paste0("Protein_", 1:n_proteins)  # Protein names


# Create corresponding labels (conditions)
labels <- c(rep("Control", 5), rep("Treatment1", 5), rep("Treatment2", 5))


sample_data <- lfq_data %>% 
  as.data.frame() 

sample_data$label <- c(rep("Control", 5), rep("Treatment1", 5), rep("Treatment2", 5))


sample_data$sample <- paste("sample_", 1:15)

sample_data <- sample_data %>% 
  select(c("sample", "label"), everything())


write.csv(sample_data, "data/sample_data.csv", row.names = FALSE)
