pan.rezultati <- read.csv("Data/Pannzer_human proteome_2.txt", header = T, 
                            sep = "\t")

pan.cc <- pan.rezultati[pan.rezultati$ontology == "CC", ]


frenki <- read.table("Data/FRENKI_CC_predikcije_propagirano.txt", header = F, 
                     sep = "\t")

samo.frenki <- as.data.frame(setdiff(frenki$V1, pan.rezultati$qpid))
names(samo.frenki) <- c("Protein")

samo.frenki.cc <- as.data.frame(setdiff(frenki$V1, pan.cc$qpid))
names(samo.frenki.cc) <- c("Protein")

frenki3 <- frenki[frenki$V1 %in% samo.frenki$Protein, ] 
frenki3 <- frenki3[frenki3$V2 == "GO:0005739", ]

frenki4 <- frenki[frenki$V1 %in% samo.frenki.cc$Protein, ] 
frenki4 <- frenki4[frenki4$V2 == "GO:0005739", ]


mapa <- read.csv("Data/CAFA_mapping.txt", sep = "\t", header = F)

frenki3 <- merge(frenki3, mapa[, c(1, 4)], by.x = "V1", 
                 by.y = "V4", all.x = T)
names(frenki3) <- c("Cafa", "go", "Score", "Protein")

frenki4 <- merge(frenki4, mapa[, c(1, 4)], by.x = "V1", 
                 by.y = "V4", all.x = T)
names(frenki4) <- c("Cafa", "go", "Score", "Protein")



impi <- read.table("Data/impi_proteini_reviewed.tab",
                   header = T, sep = "\t")
names(impi) <- c("Gen", "Protein")

provera <- as.data.frame(intersect(impi$Protein, frenki3$Protein))

sum(frenki3$Protein %in% impi$Protein)


provera4 <- as.data.frame(intersect(impi$Protein, frenki4$Protein))
 