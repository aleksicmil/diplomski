### 
### Da li su protieni iy IMPI-ja u CAFA training setu
###
###________________________________________________________________

# Ucitavanje fajlova

impi <- read.table("Data/All_mito_info.txt", sep = "\t", header = T)
cafa <- read.table("Data/trening_cafa_C_ext.txt", sep = "\t", header = F, 
                   col.names = c("Protein", "Ontology"))

impi_protein <- read.csv("Data/impi_protein.csv", header = F, 
                         col.names = "Protein")

impi_oba <- read.table("Data/impi_UniProt_to_ENSG.txt", header = T, sep = "\t")

# Pravi listu svih proteina u CAFA-i

cafa_protein <- as.data.frame(levels(cafa$Protein))
cafa_protein[, 2] <- "da"
names(cafa_protein) <- c("Protein", "CAFA")


# 


sum(impi_protein$Protein %in% levels(cafa$Protein))

match(impi_protein$Protein, cafa$Protein, nomatch = NA)



novo <- merge(impi, impi_oba, by = "Gene.ID")
novo_GO <- merge(novo, cafa, by = "Protein", all.x = T)
novo_cafa <- merge(novo, cafa_protein, by = "Protein", all.x = T)
table(novo_cafa$CAFA == "da")

table(novo_cafa$CAFA != "da")
write.table(impi$Gene.ID, "Data/impi_ensg.txt", col.names = F, row.names = F, 
            quote = F)

novo <- novo[,c(1, 54, 2:53)]
