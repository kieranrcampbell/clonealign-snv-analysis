suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(aargh))

parse_vcf <- function(input_vcf = "input.vcf",
                      cell_barcode = "cell_barcode",
                      clone = "A") {
  vcf <- readVcf(input_vcf, "hg19")
  snvs <- names(ranges(vcf))
  if(length(snvs) > 0) {
    snv_txt <- sapply(snvs, function(snv) {
      paste(cell_barcode, clone, paste0(snv, "\n"), sep = ",")
    })
    cat(paste(snv_txt, collapse = ""))
  }
}

aargh(parse_vcf)
