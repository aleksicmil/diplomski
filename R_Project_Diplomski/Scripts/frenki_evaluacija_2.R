# FRENKI EVALUACIJA ponovljeno

#_________________________________________________________________________

# UCITAVANJE DATASETOVA

# Gene ontology - mitochondrion
go <- "GO:0005739"

# Radni dataset, proteini iz IMPI-ja koji nisu u CAFA training datasetu (2827)
radni <- read.table("Data/radni_dataset.txt", sep = "\t", header = T)
radni.uni <- radni[!(duplicated(radni$Protein)), ]
radni.uni$Protein <- droplevels(radni.uni)$Protein

# Frenki predikcije. Izdvojeni proteini za koje je predvidjeno da su mito (3163)
frenki <- read.table("Data/FRENKI_predikcije_mapirano.txt", 
                     header = T, 
                     sep = "\t")
frenki.mito <- frenki[frenki$Predikcija == go, ]
frenki.mito$UniProt <- droplevels(frenki.mito)$UniProt

# Ucitavanje CAFA trening seta. Vadjenje samo mito proteina. (4675)

cafa <- read.table("Data/trening_cafa_C_ext.txt", header = F, sep = "\t")
names(cafa) <- c("UniProt", "Predikcija")
cafa.mito <- cafa[cafa$Predikcija == go, ]
names(cafa.mito) <- c("UniProt", "Predikcija")
rm(cafa)

# _________________________________________________________________________

# CISCENJE FRENKI.MITO

# Izbacivanje duplikata iz frenki.mito

dim(frenki.mito)[1] > length(levels(frenki.mito$UniProt))
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

t.pos <- as.data.frame(intersect(radni.uni$Protein, frenki.mito$UniProt)) 
names(t.pos) <- c("UniProt")
tp <- dim(t.pos)[1] # tp = 47

# False positives - proteini za koje je FRENKI rekao da su mitohondrijalni, a
#                   NE nalaze se u radnom. Ostatak frenki.mito kada se oduzme tp

f.pos <- as.data.frame(setdiff(frenki.mito$UniProt, t.pos$UniProt))
fp <- dim(f.pos)[1] # fp = 2715

# False negatives - svi oni proteini koje FRENKI nije predvideo kao mito, 
#                   a nalaze se u radnom. Ostatak radnog kada se oduzme tp

f.neg <- as.data.frame(setdiff(radni.uni$Protein, t.pos$UniProt))
fn <- dim(f.neg)[1] # fn = 2398

# Precision

prec <- tp/(tp+fp) # 0.017
rec <- tp/(tp+fn) # 0.019

f <- 2*prec*rec/(prec+rec) # 0.018

