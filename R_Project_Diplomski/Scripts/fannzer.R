# pan.rezultati <- read.csv("Data/Pannzer_human proteome_2.txt", header = T, 
#                             sep = "\t")
# 
# pan.cc <- pan.rezultati[pan.rezultati$ontology == "CC", ]
# 
# 
# frenki <- read.table("Data/FRENKI_CC_predikcije_propagirano.txt", header = F, 
#                      sep = "\t")
# 
# samo.frenki <- as.data.frame(setdiff(frenki$V1, pan.rezultati$qpid))
# names(samo.frenki) <- c("Protein")
# 
# samo.frenki.cc <- as.data.frame(setdiff(frenki$V1, pan.cc$qpid))
# names(samo.frenki.cc) <- c("Protein")
# 
# frenki3 <- frenki[frenki$V1 %in% samo.frenki$Protein, ] 
# frenki3 <- frenki3[frenki3$V2 == "GO:0005739", ]
# 
# frenki4 <- frenki[frenki$V1 %in% samo.frenki.cc$Protein, ] 
# frenki4 <- frenki4[frenki4$V2 == "GO:0005739", ]
# 
# 
# mapa <- read.csv("Data/CAFA_mapping.txt", sep = "\t", header = F)
# 
# frenki3 <- merge(frenki3, mapa[, c(1, 4)], by.x = "V1", 
#                  by.y = "V4", all.x = T)
# names(frenki3) <- c("Cafa", "go", "Score", "Protein")
# 
# frenki4 <- merge(frenki4, mapa[, c(1, 4)], by.x = "V1", 
#                  by.y = "V4", all.x = T)
# names(frenki4) <- c("Cafa", "go", "Score", "Protein")
# 
# 
# 
# impi <- read.table("Data/impi_proteini_reviewed.tab",
#                    header = T, sep = "\t")
# names(impi) <- c("Gen", "Protein")
# 
# provera <- as.data.frame(intersect(impi$Protein, frenki3$Protein))
# 
# sum(frenki3$Protein %in% impi$Protein)
# 
# 
# provera4 <- as.data.frame(intersect(impi$Protein, frenki4$Protein))


f.tp <- read.table("Data/frenki_tp.txt", header = T, sep = "\t")
f.fp <- read.table("Data/frenki_fp.txt", header = T, sep = "\t")


p.tp <- read.table("Data/pannzer_tp.txt", header = T, sep = "\t")
p.fp <- read.table("Data/pannzer_fp.txt", header = T, sep = "\t")

centar <- as.data.frame(intersect(f.tp$UniProt, p.tp$UniProt))
names(centar) <- c("UniProt")

samo.f <- f.tp[!(f.tp$UniProt %in% centar$UniProt), ]
samo.p <- as.data.frame(p.tp[!(p.tp$UniProt %in% centar$UniProt), ])
names(samo.p) <- c("UniProt")


presek.greske <- as.data.frame(intersect(f.fp$UniProt, p.fp$UniProt))

impi.info <- read.table("Data/All_mito_info.txt", header = T, sep = "\t")
impi.prot <- read.table("Data/impi_proteini_reviewed.tab", 
                        header = T, sep = "\t")

impi <- merge(impi.info, impi.prot, by.x = "Gene.ID", by.y = "Gene", 
              all.y = T)
duzine <- read.table("Data/impi_duzina.tab", sep = "\t", header = T)
impi <- merge(impi, duzine, by.x = "UniProt", by.y = "Entry", 
              all.x = T)

samo.f <- merge(samo.f, impi, by = "UniProt", all.x = T)
samo.p <- merge(samo.p, impi, by = "UniProt", all.x = T)

samo.f$Disease <- as.factor(samo.f$Disease)
samo.p$Disease <- as.factor(samo.p$Disease)

plot(samo.f$Category)
plot(samo.p$Category)

table(samo.f$Disease)
table(samo.p$Disease)


go <- read.table("Data/GO_mitochondrion.txt", header = F, sep = "\t")

samo.f <- merge(samo.f, go[, c(1, 4)], by.x = "UniProt", by.y = "V1", 
                all.x = T)
samo.p <- merge(samo.p, go[, c(1, 4)], by.x = "UniProt", by.y = "V1", 
                all.x = T)

write.table(samo.f$UniProt, "Data/frenki_tp_UniProt.txt", 
            col.names = F, row.names = F, sep = "\t", quote = F)

plot(density(samo.p$Length), main = "Broj AK u proteinu", col = "red",
     Xlab = "Duzina proteina")
lines(density(samo.f$Length), col = "blue")
lines(density(impi$Length, na.rm = T), col = "green")
legend(x = 5000, y = 0.002, legend = c("PANNZER", "FRENKI", "IMPI"), 
       fill = c("red", "blue", "green"))

plot(density(samo.p$Length), main = "Broj AK u proteinu", col = "red",
     Xlab = "Duzina proteina", xlim = c(0, 2000))
lines(density(samo.f$Length), col = "blue")
lines(density(impi$Length, na.rm = T), col = "green")
legend(x = 1500, y = 0.002, legend = c("PANNZER", "FRENKI", "IMPI"), 
       fill = c("red", "blue", "green"))

plot(density(samo.p$Length), main = "Broj AK u proteinu", col = "red",
     Xlab = "Duzina proteina", xlim = c(0, 100))
lines(density(samo.f$Length), col = "blue")
lines(density(impi$Length, na.rm = T), col = "green")
legend(x = 1500, y = 0.002, legend = c("PANNZER", "FRENKI", "IMPI"), 
       fill = c("red", "blue", "green"))

plot(density(samo.p$Length), main = "Broj AK u proteinu", col = "red",
     Xlab = "Duzina proteina", xlim = c(0, 500))
lines(density(samo.f$Length), col = "blue")
lines(density(impi$Length, na.rm = T), col = "green")
legend(x = 1500, y = 0.002, legend = c("PANNZER", "FRENKI", "IMPI"), 
       fill = c("red", "blue", "green"))


rez <- matrix(c(0.14, 0.04, 0.6, 0.31, 0.42, 0.155, 0.76, 0.64, 
                0.21, 0.063, 0.67, 0.41),
              nrow = 3, ncol = 4, byrow = T,
              dimnames = list(c("preciznost", "odziv", "F mera"), 
              c("FANNZER", "FRENKI", "PANNZER", "PANNZER_GO")))

rez <- t(rez)
rez <- rez[c(2, 3, 4, 1), ]
barplot(rez, beside = T, ylim = c(0, 1), col = c(1, 2, 3, 4))
legend("topleft", legend = c("FRENKI", "PANNZER", "PANNZER_GO", "FANNZER"), 
       cex = 0.6, fill = c(1, 2, 3, 4))
