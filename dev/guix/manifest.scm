(define packages
 (list
;; Internal development packages:
"vim"
"nss-certs"
"git"
"glibc-locales"
"htslib"
;; Snakemake control flow:
"python"
"python-pyyaml"
"snakemake"
"python-ipython"
"graphviz"
;; R-packages:
"r-minimal"
"r-genomicalignments"
"r-rmarkdown"
"r-data-table"
"r-rtracklayer"
"r-colorout"
"r-dplyr"
;; Mapping and mod-detection:
"samtools"
"bwa"
"python-ont-tombo"
"python-rpy2"
"minimap2"
"nanopolish"
"python-ont-fast5-api"
))
(specifications->manifest packages)
