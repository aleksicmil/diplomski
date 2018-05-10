impi <- read.table("Data/All_mito_info.txt", sep = "\t", header = T)
write.table(impi$Gene.ID, "Data/impi_geni.txt", sep = "\t", col.names = F, 
            quote = F, row.names = F)
impi.prot <- read.table("Data/impi_mapirano.tab", sep = "\t", header = T)
names(impi.prot) <- c("Gen", "UniProt") 
write.table(impi.prot, "Data/impi_proteini_mapirano_novo.txt", sep = "\t", 
            quote = T, col.names = T, row.names = F)
# _______________________________________

mito.go <- read.table("Data/GO_mitochondrion.txt", header = F, sep = "\t")

izbaci <- function (x){
  x <- substr(x, 11, nchar(x))
}

mito.go$V1 <- as.character(mito.go$V1)
mito.go$V1 <- izbaci(mito.go$V1)

presek <- as.data.frame(intersect(impi.prot$UniProt, mito.go$V1))

#______________________________________________________________

pan.proteom <- read.table("Data/Pannzer_human proteome_2.txt", header = T, 
                          sep = "\t")
pan.proteom <- pan.proteom[pan.proteom$ontology == "CC", ]

write.table(pan.proteom, "Data/pan_proteom_cc.txt", quote = F, 
            col.names = T, row.names = F)
