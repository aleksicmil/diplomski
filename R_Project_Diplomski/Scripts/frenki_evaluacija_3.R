# FRENKI EVALUACIJA ponovljeno

#_________________________________________________________________________

# UCITAVANJE DATASETOVA

# Gene ontology - mitochondrion
go <- "GO:0005739"

# Ucitavanje CAFA trening seta. Vadjenje samo mito proteina. (4675)

cafa <- read.table("Data/trening_cafa_C_ext.txt", header = F, sep = "\t")
names(cafa) <- c("UniProt", "Predikcija")
cafa.mito <- cafa[cafa$Predikcija == go, ]
names(cafa.mito) <- c("UniProt", "Predikcija")
rm(cafa)

# Radni dataset, proteini iz IMPI-ja koji nisu u CAFA 
# training datasetu (1547)
radni <- read.table("Data/impi_proteini_reviewed.tab", sep = "\t", 
                    header = T)
names(radni) <- c("Gen", "UniProt")
# nema duplikata

# CIscenje radnog od kafe (710 novi radni impi minus cafa)
radni.cafa <- as.data.frame(intersect(radni$UniProt, cafa.mito$UniProt))
names(radni.cafa) <- c("UniProt")
radni.cisto <- radni[!(radni$UniProt %in% radni.cafa$UniProt), ]

# Frenki predikcije. Izdvojeni proteini za koje je predvidjeno 
# da su mito (3163, a kada se ocisti kao dole u kodu ostaje 24762 predikcija)

frenki <- read.table("Data/FRENKI_predikcije_mapirano.txt", 
                     header = T, 
                     sep = "\t")
frenki.mito <- frenki[frenki$Predikcija == go, ]
frenki.mito$UniProt <- droplevels(frenki.mito)$UniProt
write.table(frenki.mito, "Data/frenki.mito.txt", quote = F, col.names = T, 
            row.names = F, sep = "\t")
# _________________________________________________________________________

# CISCENJE FRENKI.MITO

# Izbacivanje duplikata iz frenki.mito

dim(frenki.mito)[1] - length(levels(frenki.mito$UniProt))
frenki.mito <- frenki.mito[!duplicated(frenki.mito$UniProt), ]

# Izbacivanje NA za UniProt iz frenki.mito

frenki.mito <- frenki.mito[!(is.na(frenki.mito$UniProt)), ]

# Izbacivanje proteina koji su bili u CAFA trening setu iz frenki.mito

presek <- as.data.frame(intersect(frenki.mito$UniProt, cafa.mito$UniProt))
names(presek) <- c("UniProt")

frenki.mito <- frenki.mito[!(frenki.mito$UniProt %in% cafa.mito$UniProt), ]

#__________________________________________________________________________

# EVALUACIJA

# True positives - broj proteina za koje je FRENKI predvideo da su mito, 
#                  a nalaze se u radnom datasetu

t.pos <- as.data.frame(intersect(radni.cisto$UniProt, frenki.mito$UniProt)) 
names(t.pos) <- c("UniProt")
tp <- dim(t.pos)[1] # tp = 110

# False positives - proteini za koje je FRENKI rekao da su mitohondrijalni, a
#                   NE nalaze se u radnom. Ostatak frenki.mito kada se oduzme tp

f.pos <- as.data.frame(setdiff(frenki.mito$UniProt, t.pos$UniProt))
names(f.pos) <- c("UniProt")
fp <- dim(f.pos)[1] # fp = 2652

# False negatives - svi oni proteini koje FRENKI nije predvideo kao mito, 
#                   a nalaze se u radnom. Ostatak radnog kada se oduzme tp

f.neg <- as.data.frame(setdiff(radni.cisto$UniProt, t.pos$UniProt))
names(f.neg) <- c("UniProt")
fn <- dim(f.neg)[1] # fn = 600

# Precision

prec <- tp/(tp+fp) # 0.040
rec <- tp/(tp+fn) # 0.155

f <- 2*prec*rec/(prec+rec) # 0.063

# ___________________________________________

# Pisanje resenja

write.table(f.tp, "Data/frenki_tp.txt", row.names = F, col.names = T, 
            quote = F, sep =  "\t")
write.table(f.fp, "Data/frenki_fp.txt", row.names = F, col.names = T, 
            quote = F, sep =  "\t")


# Crtanje dijagrama

draw.pairwise.venn(710, 2762, 110, 
                   category = c("FRENKI predikije", "IMPI - CAFA"), 
                   scaled = T, inverted = T)

# Priprema fajlova za poredjenja algoritama

# Pravi datasetove za tp, fp i fn

f.tp <- radni.cisto[radni.cisto$UniProt %in% t.pos$UniProt, ]
f.tp[, 3] <- "tp"
names(f.tp)[3] <- "Guess"

f.fn <- radni.cisto[radni.cisto$UniProt %in% f.neg$UniProt, ]
f.fn[, 3] <- "fn"
names(f.fn)[3] <- "Guess"

f.res <- rbind(f.tp, f.fn)
f.res[ ,4] <- "FRENKI"
names(f.res)[4] <- c("Alg") 

write.table(f.res, "Data/frenki_resenja_tpfn.txt", quote = F, row.names = F, 
            col.names = T, sep = "\t")

# False positivi se ne nalaze u radnom datasetu, ne mogu se spojiti sa tp i fn

f.fp <- frenki.mito[frenki.mito$UniProt %in% f.pos$UniProt, ]
f.fp[, 3] <- "fp"
names(f.fp)[7] <- "Guess"

write.table(f.fp, "Data/frenki_resenja_fp.txt", sep = "\t", quote = F, 
            col.names = T, row.names = F)
