###
### Pravljenje dataseta za dalji rad i mapiranje predikcija
###
###________________________________________________________


# impi sa mitominera
impi <- read.table("Data/impi_mitominer.tsv", sep = "\t", header = F, 
                   col.names = c("Protein", "Gene.ID", "Gene", "Status"))
# svi geni u impiju
impi_genes <- as.data.frame(levels(impi$Gene.ID))

# cafa training set
cafa <- read.table("Data/trening_cafa_C_ext.txt", header = F, sep = "\t", 
                   col.names = c("Protein", "Ontology"))

# lista svih proteina u cafi posto se u cafa ponavljaju zbog ontologija
cafa_proteini <- as.data.frame(levels(cafa$Protein))
cafa_proteini[, 2] <- cafa_proteini[, 2] <- "Da"
names(cafa_proteini) <- c("Protein", "CAFA")

# listi impi proteina dodaje informaciju da li su u cafa training datasetu
impi_u_cafi <- merge(impi, cafa_proteini, by = "Protein", all.x = T)
table(impi_u_cafi$CAFA == "Da")
write.table(impi_u_cafi, "Data/impi_info_u_cafi.txt", sep = "\t", quote = F, 
            col.names = T, row.names = F)

# pravi listu proteina iy impija koji su u cafa training daatsetu
u_cafi <- subset(impi_u_cafi, impi_u_cafi$CAFA == "Da")
write.table(u_cafi, "Data/samo_impi_u_cafi.txt", sep = "\t", quote = F, 
            col.names = T, row.names = T)

# pravi radni data set - lista mitohondrijalnih proteina iz impija koji nisu u 
# CAFA training setu

radni <- subset(impi_u_cafi, is.na(impi_u_cafi$CAFA == "Ne"))
write.table(radni, "Data/radni_dataset.txt", sep = "\t", quote = F, 
            col.names = T, row.names = T)


