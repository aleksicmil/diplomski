
# Ucitavanje PANNZer rezultata i parsovanje imena

pan <- read.csv("Data/propagated_list_pannzer_mapped.txt", sep = "\t", 
                header = T)

# _________________________________________________

# Ekstrahovanje mito predikcija iz PANNZER rezultata

pan.mito <- pan[pan$go == 5739, ]
pan.mito$UniProt <- as.factor(pan.mito$UniProt)
pan.mito$UniProt <- droplevels(pan.mito)$UniProt

#__________________________________________________

# Evaluacija

radni <- read.table("Data/impi_proteini_reviewed.tab", header = T, sep = "\t")
radni <- radni[!duplicated(radni$UniProt), ]

# True positives - proteini koji su i u PANNZERu i u radnom

t.pos <- as.data.frame(intersect(radni$UniProt, pan.mito$UniProt))
tp <- dim(t.pos)[1] #1169
names(t.pos) <- c("UniProt")

# False positives - proteini koji nisu u radnom a PANNZER je 
# predvideo da su mito

f.pos <- as.data.frame(setdiff(pan.mito$UniProt, radni$UniProt))
fp <- dim(f.pos)[1] #759
names(f.pos) <- c("UniProt")

# False negatives - sve sto PANNZER nije predvido kao mito a u radnom su

f.neg <- as.data.frame(setdiff(radni$UniProt, pan.mito$UniProt))
fn <- dim(f.neg)[1] #378
names(f.neg) <- c("uniProt")
# Precision, recall, i f score

prec <- tp/(tp+fp) # 0.6
rec <- tp/(tp+fn) # 0.76
f <- 2*prec*rec/(prec+rec) #0.67

# ____________________________________________________________________________

# Koliko je precizan premo onima koji u GO nisu obelezeni terminom GO:0005739

go <- read.table("Data/GO_mitochondrion.txt", header = F, sep = "\t")

radni.go <- radni[!(radni$UniProt %in% go$V1), ]

t.pos.go <- as.data.frame(intersect(radni.go$UniProt, pan.mito$UniProt))
tp.go <- dim(t.pos.go)[1] #593
names(t.pos.go) <- c("UniProt")

# False positives - proteini koji nisu u radnom a PANNZER je 
# predvideo da su mito

f.pos.go <- as.data.frame(setdiff(pan.mito$UniProt, radni.go$UniProt))
fp.go <- dim(f.pos.go)[1] #1335
names(f.pos.go) <- c("UniProt")

# False negatives - sve sto PANNZER nije predvido kao mito a u radnom su

f.neg.go <- as.data.frame(setdiff(radni.go$UniProt, pan.mito$UniProt))
fn.go <- dim(f.neg.go)[1] #337
names(f.neg.go) <- c("UniProt")
# Precision, recall, i f score

prec.go <- tp.go/(tp.go+fp.go) # 0.31
rec.go <- tp.go/(tp.go+fn.go) # 0.64
f.go <- 2*prec.go*rec.go/(prec.go+rec.go) #0.41

# Pravljenjej fajlau kome su sve tačne predikcije pannzera

write.table(t.pos.go, "Data/pannzer_tp.txt", row.names = F, col.names = T, 
            sep = "\t")
write.table(f.pos.go, "Data/pannzer_fp.txt", row.names = F, col.names = T, 
            sep = "\t")


# # Priprema za poredjenje algoritama
# 
# # Pravi datasetove za tp, fp i fn
# 
# p.tp <- radni[radni$Protein %in% t.pos$Protein, ]
# p.tp[ ,6] <- "tp"
# names(p.tp)[6] <- c("Guess")
# 
# p.fn <- radni[radni$Protein %in% f.neg$Protein, ]
# p.fn[ ,6] <- "fn"
# names(p.fn)[6] <- c("Guess")
# 
# p.res1 <- rbind(p.tp, p.fn)
# p.res1 <- p.res1[ ,-5]
# p.res1[, 6] <- "PANNZER"
# names(p.res1)[6] <- c("Alg")
# write.table(p.res1, "Data/pannzer_resenja_tpfn.txt", quote = F, sep = "\t", 
#             col.names = T, row.names = F)
# 
# # False positivi ne postoje u radno datasetu i zato ne mogu da ih spojim sa
# # onim resenjima
# 
# p.fp <- pan.mito[pan.mito$Protein %in% f.pos$Protein, ]
# p.fp[ ,9] <- "fp"
# names(p.fp)[9] <- c("Guess")
# 
# write.table(p.fp, "Data/pannzer_resenja_fp.txt", sep = "\t", quote = F, 
#             col.names = T, row.names = F)

