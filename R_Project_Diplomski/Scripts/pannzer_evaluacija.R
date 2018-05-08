
# Ucitavanje PANNZer rezultata i parsovanje imena

pan <- read.table("Data/pannzer_go_rezultati.txt", sep = "\t", header = T)
pan$qpid <- as.character(panz$qpid)

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
tp <- dim(t.pos)[1]
names(t.pos) <- c("Protein")

# False positives - proteini koji nisu u radnom a PANNZER je predvideo da su mito

f.pos <- as.data.frame(setdiff(pan.mito$Protein, radni$Protein))
fp <- dim (f.pos)[1]

# False negatives - sve sto PANNZER nije predvido kao mito a u radnom su

f.neg <- as.data.frame(setdiff(radni$Protein, pan.mito$Protein))
fn <- dim(f.neg)[1]

# Precision, recall, i f score

prec <- tp/(tp+fp)
rec <- tp/(tp+fn)
f <- 2*prec*rec/(prec+rec)

pog <- radni[radni$Protein %in% t.pos$Protein, ]

barplot(table(pog$Status))
