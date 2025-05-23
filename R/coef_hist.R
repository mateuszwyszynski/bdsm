#' Graphs of the distribution of the coefficients over the model space
#'
#' This function draws graphs of the distribution (in the form of histogram or kernel density) of the coefficients for all the considered regressors over the part of the model space that includes this regressors (half of the model space).
#'
#' @name coef_hist
#'
#' @param bma_list bma object (the result of the bma function)
#' @param BW Parameter indicating what method should be chosen to find bin widths for the histograms: \cr
#' 1) "FD" Freedman-Diaconis method \cr
#' 2) "SC" Scott method \cr
#' 3) "vec" user specified bin widths provided through a vector (parameter: binW)
#' @param binW A vector with bin widths to be used to construct histograms for the regressors. The vector must be of the size equal to total number of regressors. The vector with bin widths is used only if parameter BW="vec".
#' @param BN Parameter taking the values (default: BN = 0): \cr
#' 1 - the histogram will be build based on the number of bins specified by the user through parameter num. If BN=1, the function ignores parameters BW. \cr
#' 0 - the histogram will be build in line with parameter BW
#' @param num A vector with the numbers of bins used to be used to construct histograms for the regressors. The vector must be of the size equal to total number of regressors. The vector with bin widths is used only if parameter BN=1.
#' @param kernel A parameter taking the values (default: kernel = 0):\cr
#' 1 - the function will build graphs using kernel density for the distribution of coefficients (with kernel=1, the function ignores parameters BW and BN) \cr
#' 0 - the function will build regular histogram density for the distribution of coefficients
#'
#' @return A list with the graphs of the distribution of coefficients for all the considered regressors.
#'
#' @export
#'
#' @examples
#' \donttest{
#' library(magrittr)
#'
#' data_prepared <- bdsm::economic_growth[, 1:6] %>%
#'   bdsm::feature_standardization(
#'     excluded_cols = c(country, year, gdp)
#'   ) %>%
#'   bdsm::feature_standardization(
#'     group_by_col  = year,
#'     excluded_cols = country,
#'     scale         = FALSE
#'   )
#'
#' bma_results <- bma(
#'   model_space = bdsm::small_model_space,
#'   df          = data_prepared,
#'   round       = 3,
#'   dilution    = 0
#' )
#'
#' coef_plots <- coef_hist(bma_results, kernel = 1)
#' }

utils::globalVariables(".data")

coef_hist <- function(bma_list, BW = "FD", binW = NULL, BN = 0, num = NULL, kernel = 0){

x_names <- bma_list[[3]] # names of variables
K <- bma_list[[4]] + 1 # number of variables
alpha <- bma_list[[13]]
betas <- bma_list[[14]]

# Adding colnames and changing to dataframe
colnames(alpha) <- x_names[1]
colnames(betas) <- x_names[-1]
alpha <- as.data.frame(alpha)
betas <- as.data.frame(betas)

histPlots<-list() # Opening a list  for the histogram plots

# CONDITION for using kernel or regular histogram
if (kernel==1){
  histPlots[[1]]<-invisible(ggplot2::ggplot(alpha, ggplot2::aes(x = .data[[x_names[1]]])) +
                              ggplot2::geom_density(fill = "skyblue", alpha = 0.7) +
                              ggplot2::labs(
                                title = paste("Distribiution of", x_names[1], "coefficients"),
                                x = paste0("Coefficients on ",x_names[1]),
                                y = "Frequency") +
                              ggplot2::theme_minimal(base_size = 12) +
                              ggplot2::theme(
                                plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                axis.title = ggplot2::element_text(face = "bold")))
  names(histPlots)[[1]] <-x_names[[1]]
  for (i in 2:K){
    histPlots[[i]]<-invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                ggplot2::geom_density(fill = "skyblue", alpha = 0.7) +
                                ggplot2::labs(
                                  title = paste("Distribiution of", x_names[i], "coefficients"),
                                  x = paste0("Coefficients on ",x_names[i]),
                                  y = "Frequency") +
                                ggplot2::theme_minimal(base_size = 12) +
                                ggplot2::theme(
                                  plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                  axis.title = ggplot2::element_text(face = "bold")))
    names(histPlots)[[i]] <-x_names[[i]]
  }

}else{# REGULAR HISTOGRAM BELOW

  # CONDITION for graphs plotted with binwidth:
  if (BN==0){### Rules for bin width
    # 1) Freedman-Diaconis (FD)
    if (BW=="FD"){BW<-(stats::IQR(alpha[,1])*2)/sqrt(length(alpha[,1]))}
    # 2) Scott (SC)
    if (BW=="SC"){BW<-(stats::sd(alpha[,1])*3.5)/(length(alpha[,1])^(1/3))}
    # 3) Binwidth sizes
    if(BW=="vec"){
      if (is.null(binW)){stop("Please provide a vector with bin width sizes through parameter binW")}
      if (length(binW)!=K){stop("binW is missspecified: binW should have K (number of regressors +1) elements")}
      BW<-binW[1]
    }
    histPlots[[1]] <- invisible(ggplot2::ggplot(alpha, ggplot2::aes(x = .data[[x_names[1]]])) +
                                  ggplot2::geom_histogram(binwidth = BW, fill = "skyblue", color = "white", alpha = 0.8) +
                                  ggplot2::labs(
                                    title = paste("Distribiution of", x_names[1], "coefficients"),
                                    x = paste0("Coefficients on ",x_names[1]),
                                    y = "Frequency") +
                                  ggplot2::theme_minimal(base_size = 12) +
                                  ggplot2::theme(
                                    plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                    axis.title = ggplot2::element_text(face = "bold")))
    names(histPlots)[[1]] <- x_names[[1]]
    for (i in 2:K){ # at this LOOP we go through all the regressors
      # 1) Freedman-Diaconis (FD)
      if (BW=="FD"){BW<-(stats::IQR(betas[,i])*2)/sqrt(length(betas[,i]))}
      # 2) Scott (SC)
      if (BW=="SC"){BW<-(stats::sd(betas[,i])*3.5)/(length(betas[,i])^(1/3))}
      # 3) Binwidth sizes
      if(BW=="vec"){
        if (is.null(binW)){stop("Please provide a vector with bin width sizes through parameter binW")}
        if (length(binW)!=K){stop("binW is missspecified: binW should have K elements")}
        BW<-binW[i]
      }
      histPlots[[i]] <- invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                    ggplot2::geom_histogram(binwidth=BW, fill = "skyblue", color = "white", alpha = 0.8) +
                                    ggplot2::labs(
                                      title = paste("Distribiution of", x_names[i], "coefficients"),
                                      x = paste0("Coefficients on ",x_names[i]),
                                      y = "Frequency") +
                                    ggplot2::theme_minimal(base_size = 12) +
                                    ggplot2::theme(
                                      plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                      axis.title = ggplot2::element_text(face = "bold")))
      names(histPlots)[[i]] <- x_names[[i]]
    }
  }

  # CONDITION for graphs plotted with bins - through setting the number of bins:
  if (BN==1){### Rules for bin width
    if (is.null(num)){stop("Please provide a vector with number of bins through parameter num")}
    if (length(num)!=K){stop("num is missspecified: num should have K elements")}
    histPlots[[1]]<-invisible(ggplot2::ggplot(alpha, ggplot2::aes(x = .data[[x_names[1]]])) +
                                ggplot2::geom_histogram(bins=num[1], fill = "skyblue", color = "white", alpha = 0.8) +
                                ggplot2::labs(
                                  title = paste("Distribiution of", x_names[1], "coefficients"),
                                  x = paste0("Coefficients on ",x_names[1]),
                                  y = "Frequency") +
                                ggplot2::theme_minimal(base_size = 12) +
                                ggplot2::theme(
                                  plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                  axis.title = ggplot2::element_text(face = "bold")))
    names(histPlots)[[1]] <-x_names[[1]]
    for (i in 2:K){
      histPlots[[i]]<-invisible(ggplot2::ggplot(betas, ggplot2::aes(x = .data[[x_names[i]]])) +
                                  ggplot2::geom_histogram(bins=num[i], fill = "skyblue", color = "white", alpha = 0.8) +
                                  ggplot2::labs(
                                    title = paste("Distribiution of", x_names[i], "coefficients"),
                                    x = paste0("Coefficients on ",x_names[i]),
                                    y = "Frequency") +
                                  ggplot2::theme_minimal(base_size = 12) +
                                  ggplot2::theme(
                                    plot.title = ggplot2::element_text(size = 12, hjust = 0.5, face = "bold"),
                                    axis.title = ggplot2::element_text(face = "bold")))
      names(histPlots)[[i]] <-x_names[[i]]
    }
  }
}
return(histPlots)
}
