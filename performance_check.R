# CASE 1: SEM_likelihood when transformation of data and params is needed.
microbenchmark::microbenchmark({
  set.seed(1)
  SEM_likelihood(
    0.5,
    generate_test_feature_standard_data(),
    times, entities, dep_var
  )
})
# Unit: milliseconds
# min      lq       mean.    median   uq       max        neval
# 25.80372 26.05294 27.14582 26.22264 27.25252 36.02822   100

# (data and params in the form ready for SEM_likelihood)
set.seed(1)
data <- generate_test_feature_standard_data()
data <- data %>%
  matrices_from_df(timestamp_col = times,
                   entity_col = entities,
                   dep_var_col = dep_var)
params <- SEM_params_to_list(0.5, data)
# CASE 2: SEM_likelihood when transformation of data and params is NOT needed,
# but CHECKED.
microbenchmark::microbenchmark({
  set.seed(1)
  SEM_likelihood(
    params,
    data,
    times, entities, dep_var
  )
}, unit = "milliseconds")
# Unit: milliseconds
# min        lq        mean       median    uq        max        neval
# 0.064206   0.065108  0.06850403 0.0656    0.0666865 0.200326   100

# CASE 3: SEM_likelihood when transformation of data and params is NOT needed,
# and NOT CHECKED.
microbenchmark::microbenchmark({
  set.seed(1)
  SEM_likelihood_bma(
    params,
    data,
    times, entities, dep_var
  )
}, unit = "milliseconds")
# Unit: milliseconds
# min        lq        mean       median    uq        max        neval
# 0.057564   0.0583635 0.06094978 0.0591015 0.0597165 0.182983   100
