library(adegenet)
library(recipes)
library(assignPOP)

ref <- read.Structure("06_newsamples/reference_5pop.stru")
new <- read.Structure("06_newsamples/unknown.stru")

assign.X(x1=ref,
         x2=new, dir="06_newsamples/assignment_5pop/",
         model="svm")
