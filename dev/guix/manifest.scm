(define packages 
 (list
 ;; Internal development:
"vim"
"git"
"glibc-locales"
;; Snakemake control flow:
"python"
"snakemake"
"python-ipython"
"python-pyyaml"
"graphviz"
"htslib"
;; R-packages:
"r-minimal"
"r-genomicalignments"
"r-rmarkdown"
"r-data-table"
"r-rtracklayer"
"r-dplyr"
;; Mapping and mod-detection:
"samtools"
"python-ont-tombo"
"minimap2"
"nanopolish"
))
(specifications->manifest packages)
