# Dispatcher wrappers for pathway analyses (supports 1..3 methods)
# - DEG functional dispatcher run_DEG_functional_analysis
# - "ClusterProfiler" used as GO
# - if method is NULL or length 0 returns invisible(NULL)

run_DEG_functional_analysis <- function(genes,
                                        method = NULL,
                                        ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "Metascape", "gProfiler2"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("ClusterProfiler_DEG_functional_analysis", mode = "function")) {
               ClusterProfiler_DEG_functional_analysis(genes = genes, ...)
             } else stop("ClusterProfiler_DEG_functional_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_DEG_functional_analysis", mode = "function")) {
               Metascape_DEG_functional_analysis(genes = genes, ...)
             } else stop("Metascape_DEG_functional_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_DEG_functional_analysis", mode = "function")) {
               gProfiler2_DEG_functional_analysis(genes = genes, ...)
             } else stop("gProfiler2_DEG_functional_analysis not found")
           },
           stop("Unknown DEG functional method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}

run_overrepresentation_analysis <- function(genes,
                                           method = NULL,
                                           ...) {
  if (is.null(method) || length(method) == 0) return(invisible(NULL))

  method <- match.arg(method, choices = c("ClusterProfiler", "gProfiler2", "Metascape"), several.ok = TRUE)
  if (length(method) > 3) stop("up to 3 methods allowed")

  res <- lapply(method, function(m) {
    switch(m,
           ClusterProfiler = {
             if (exists("ClusterProfiler_overrepresentation_analysis", mode = "function")) {
               ClusterProfiler_overrepresentation_analysis(genes = genes, ...)
             } else stop("ClusterProfiler_overrepresentation_analysis not found")
           },
           gProfiler2 = {
             if (exists("gProfiler2_overrepresentation_analysis", mode = "function")) {
               gProfiler2_overrepresentation_analysis(genes = genes, ...)
             } else stop("gProfiler2_overrepresentation_analysis not found")
           },
           Metascape = {
             if (exists("Metascape_overrepresentation_analysis", mode = "function")) {
               Metascape_overrepresentation_analysis(genes = genes, ...)
             } else stop("Metascape_overrepresentation_analysis not found")
           },
           stop("Unknown overrepresentation method: ", m))
  })
  names(res) <- method
  if (length(res) == 1) res[[1]] else res
}
