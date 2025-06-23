print_hashes <- function(vars) {
  hashes <- sapply(vars, function(x) substr(xxhashlite::xxhash(x), 1, 4))
  cat(paste(names(hashes), hashes, collapse = " "), "\n")
}

print_debug <- function(fn_name, hits, misses) {
  all <- hits + misses
  if (all %% 1000 == 0) {
    hits_pct <- round(100 * hits / all)
    cat(glue::glue("{fn_name}: {hits} hits out of {all} ({hits_pct}%)"), "\n")
  }
}

dumb_memo <- function(fn, fn_name = NULL, debug = FALSE) {
  hits <- 0
  misses <- 0
  last_args <- NULL
  last_result <- NULL
  function(...) {
    args <- list(...)
    if (is.null(last_args) || !identical(args, last_args)) {
      last_args <<- args
      last_result <<- fn(...)
      misses <<- misses + 1
    } else {
      hits <<- hits + 1
    }
    if (debug) {
      print_debug(fn_name, hits, misses)
    }
    last_result
  }
}

memoise <- function(fn, debug = FALSE) {
  fn_name <- substr(xxhashlite::xxhash(substitute(fn)), 1, 4)
  memo_fn <- dumb_memo(fn, fn_name, debug = debug)
  call_fn <- function(...) {
    mc <- match.call()
    mc[[1]] <- memo_fn
    eval(mc, parent.frame())
  }
  formals(call_fn) <- formals(args(fn))
  call_fn
}
