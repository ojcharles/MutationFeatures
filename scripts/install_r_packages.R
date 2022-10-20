# this simply installs a set of packages we need

cran_packages = c("Peptides", "FSelector", "protr", "bio3d")
for( p  in cran_packages){
    install.packages(p)
}


bioc_packages = c("Biostrings")
for( p  in bioc_packages){
    BiocManager::install(p)
}


