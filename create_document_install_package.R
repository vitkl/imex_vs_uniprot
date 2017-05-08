# install devtools and roxygen2
install.packages("devtools", "roxygen")
library(devtools)
devtools::install_github("klutometis/roxygen")
library(roxygen2)

# create a package
setwd("parent_directory")
create("cats")

# document package
library(roxygen2)
setwd("/Users/vitalii/Desktop/clusterProfiler")
setwd("/Users/vitalii/Desktop/DOSE")
document()

#install package
library(devtools)
setwd("..")
install("clusterProfiler")
install("DOSE")

library(clusterProfiler)
