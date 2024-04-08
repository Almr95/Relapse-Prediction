extreme_quant_norm <- function(x,quant=0.01) {
  (x - quantile(x,quant)) / (quantile(x,1-quant) - quantile(x,quant))
}

