
# Ucitavanje PANNZer rezultata i parsovanje imena

pan <- read.table("Data/pannzer_go_rezultati.txt", sep = "\t", header = T)
pan$qpid <- as.character(pan$qpid)

ime <- function(x) {
  m1 <- gregexpr("|", x, fixed = T)
  substr(x, m1[[1]][1]+1, m1[[1]][2]-1)
}

for(i in 1:length(pan$qpid)){
  pan$Protein[i] <- ime(pan$qpid[i])
}
pan <- pan[, c(1, 8, 2:7)]

go <- 5739

# _________________________________________________

# Ekstrahovanje mito predikcija iz PANNZER rezultata

pan.mito <- pan[pan$goid == go, ]
pan.mito$Protein <- as.factor(pan.mito$Protein)
pan.mito$Protein <- droplevels(pan.mito)$Protein

#__________________________________________________

# Evaluacija

radni <- read.table("Data/radni_dataset.txt", header = T, sep = "\t")
radni <- radni[!duplicated(radni$Protein), ]

# True positives - proteini koji su i u PANNZERu i u radnom

t.pos <- as.data.frame(intersect(radni$Protein, pan.mito$Protein))
tp <- dim(t.pos)[1] #230
names(t.pos) <- c("Protein")

# False positives - proteini koji nisu u radnom a PANNZER je predvideo da su mito

f.pos <- as.data.frame(setdiff(pan.mito$Protein, radni$Protein))
fp <- dim (f.pos)[1] #1
names(f.pos) <- c("Protein")

# False negatives - sve sto PANNZER nije predvido kao mito a u radnom su

f.neg <- as.data.frame(setdiff(radni$Protein, pan.mito$Protein))
fn <- dim(f.neg)[1] #2215
names(f.neg) <- c("Protein")
# Precision, recall, i f score

prec <- tp/(tp+fp) # 0.99
rec <- tp/(tp+fn) # 0.094
f <- 2*prec*rec/(prec+rec) #0.17

# ____________________________________________________________________________

# Priprema za poredjenje algoritama

# Pravi datasetove za tp, fp i fn

p.tp <- radni[radni$Protein %in% t.pos$Protein, ]
p.tp[ ,6] <- "tp"
names(p.tp)[6] <- c("Guess")

p.fn <- radni[radni$Protein %in% f.neg$Protein, ]
p.fn[ ,6] <- "fn"
names(p.fn)[6] <- c("Guess")

p.res1 <- rbind(p.tp, p.fn)
p.res1 <- p.res1[ ,-5]
p.res1[, 6] <- "PANNZER"
names(p.res1)[6] <- c("Alg")
write.table(p.res1, "Data/pannzer_resenja_tpfn.txt", quote = F, sep = "\t", 
            col.names = T, row.names = F)

# False positivi ne postoje u radno datasetu i zato ne mogu da ih spojim sa
# onim resenjima

p.fp <- pan.mito[pan.mito$Protein %in% f.pos$Protein, ]
p.fp[ ,9] <- "fp"
names(p.fp)[9] <- c("Guess")

write.table(p.fp, "Data/pannzer_resenja_fp.txt", sep = "\t", quote = F, 
            col.names = T, row.names = F)

