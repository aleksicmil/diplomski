
pan <- read.table("Data/pannzer_go_rezultati.txt", sep = "\t", header = T)
panz <- pan
panz$qpid <- as.character(panz$qpid)

pan.cc <- subset(pan, pan$ontology == "CC")

pan.mito <- pan.cc[pan.cc$goid == "5739", ]

for (i in c(1: length(panz$qpid))){
  a = nchar(panz[i, 1])
  panz[i, 1] = substr(panz[i, 1], 11, a)
}
