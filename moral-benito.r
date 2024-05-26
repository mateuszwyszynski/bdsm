# The following packages are not installed automatically by renv,
# because the panels package does not depend on them.
# You should install them before running this script.
library(tidyverse)
library(readxl)
devtools::load_all()

set.seed(20)

data_prepared <- panels::economic_growth[,1:7] %>%
  feature_standardization(timestamp_col = year, entity_col = country) %>%
  feature_standardization(timestamp_col = year, entity_col = country,
                          cross_sectional = TRUE, scale = FALSE)

regressors <- regressor_names(data_prepared, year, country, gdp)

begin<-Sys.time()

bma_result <- bma_summary(df = data_prepared, dep_var_col = gdp,
                          timestamp_col = year, entity_col = country,
                          model_space = economic_growth_ms,
                          projection_matrix_const = TRUE)

print(paste("Computation Time:", Sys.time()-begin))

parameters_summary(regressors = regressors,
                   bet = bma_result$bet, pvarh = bma_result$pvarh,
                   pvarr = bma_result$pvarr, fy = bma_result$fy,
                   fyt = bma_result$fyt, ppmsize = bma_result$ppmsize,
                   cout = bma_result$cout, nts = bma_result$nts,
                   pts = bma_result$pts, variables_n = bma_result$variables_n)
