
# Lista _gena_ koji se nalaye u IMPIju.

impi <- read.table("Data/All_mito_info.txt", sep = "\t", header = T)
write.table(impi$Gene.ID, "Data/impi_geni.txt", sep = "\t", col.names = F, 
            quote = F, row.names = F) #ens anotacije

# Ucitam fajl za koji sada nisam sigurna odakle je ali sadrzi mapirane
# ens anotacije za IMPI prevedene u uniprot

impi.prot <- read.table("Data/impi_mapirano.tab", sep = "\t", header = T)
names(impi.prot) <- c("Gen", "UniProt") 
write.table(impi.prot, "Data/impi_proteini_mapirano_novo.txt", sep = "\t", 
            quote = T, col.names = T, row.names = F)
# _______________________________________

# Pretpostavljam da je ovo skinuti sve proteini sa Go koji imaju 
# GO mitochondrion

mito.go <- read.table("Data/GO_mitochondrion.txt", header = F, sep = "\t")

# Pravi funkciju koja od generickog GO imena pravi UniProt anotaciju
izbaci <- function (x){
  x <- substr(x, 11, nchar(x))
}

mito.go$V1 <- as.character(mito.go$V1)
mito.go$V1 <- izbaci(mito.go$V1)

# Pravi listu proteina koji su u impiju i u go kao mitohondrijalni

impi.reviewed <- read.table("Data/impi_proteini_reviewed.tab", sep ="\t",
                            header = T)
names(impi.reviewed) <- c("Gene", "UniProt")
presek <- as.data.frame(intersect(impi.reviewed$UniProt, mito.go$V1))

#______________________________________________________________

# Ovo bi trebalo daje reyultat koji se dobije kad se kroy panyer propusti 
# ceo ljudski proteom

pan.proteom <- read.table("Data/Pannzer_human proteome_2.txt", header = T, 
                          sep = "\t")
pan.proteom <- pan.proteom[pan.proteom$ontology == "CC", ]

write.table(pan.proteom, "Data/pan_proteom_cc.txt", quote = F, 
            col.names = T, row.names = F)
