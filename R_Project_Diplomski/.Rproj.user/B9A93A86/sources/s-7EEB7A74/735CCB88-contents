
###_________________________________________________________
### Provera za FRENKIja - Precision and recall
###_________________________________________________________

## Gene Ontologija mitochondrion
go <- "GO:0005739"

## Proteini iz IMPIja koji _nisu_ u CAFA training setu
radni <- read.table("Data/radni_dataset.txt", sep = "\t", header = T)
#table(radni$Status)

# ## Predvidjeni i poznati proteini u IMPIJU (svi)
# impi <- read.table("Data/impi_mitominer.tsv", header = F, sep = "\t", quote = "")
# table(impi$V4)

## FRENKI predikcije, propagirano do najstarije ontologije, mapirano na UniProt
frenki <- read.table("Data/FRENKI_predikcije_mapirano.txt", sep = "\t", 
                     header = T)

## Lista svih proteina koji postoje u FRENKI predikcijama
frenki_proteini<- as.data.frame(levels(frenki$UniProt))
names(frenki_proteini) <- c("UniProt")

## Da li se u FRENKI listi koju sam upravu ucitala nalaze proteini iz 
## training dataseta?

cafa <- read.table("Data/trening_cafa_C_ext.txt", sep = "\t", header = F)
names(cafa) <- c("UniProt", "GO")
cafa_proteini <- as.data.frame(levels(cafa[, 1]))
names(cafa_proteini) <- c("UniProt")

presek <- merge(frenki_proteini, cafa_proteini, by = "UniProt")
        #presek nije nula, tako da ima, moram da ih izbacim

## Iz FRENKIja brise one koji su u CAFA training setu
frenki_pravo <- frenki[!(frenki$UniProt %in% cafa$UniProt),]
frenki_pravo$UniProt <- as.factor(as.character(frenki_pravo$UniProt))
frenki_pravo_proteini <- as.data.frame(levels(frenki_pravo$UniProt))
  


# radni_u_frenkiju <- merge(frenki, radni, by.x = "UniProt", by.y = "Protein")
# 
# 
# frenki_mito <- frenki[frenki$Predikcija == go,]
# frenki_pravi_mito <- frenki_pravo[frenki_pravo$Predikcija == go,]

## True positives <- i FRENKI i IMPI tvrde da je protein mitohondrijalni
tp <- merge(radni, frenki, by.x = "Protein", by.y = "UniProt")
tp <- tp[tp$Predikcija == go,]

tp1 <- merge(radni, frenki_pravo, by.x = "Protein", by.y = "UniProt")
tp1 <- tp1[tp1$Predikcija == go,]

## False positives <- u FRENKIju je mito a nema ga u radnom
fp <- frenki[frenki$Predikcija == go & !(frenki$UniProt %in% radni$Protein),]

fp1 <- frenki_pravo[frenki_pravo$Predikcija == go & 
                  !(frenki_pravo$UniProt %in% radni$Protein),]


## False negatives <- ima ga u IMPIju a za FRENKIja nije mitohondrijalni
## (ali postoji u FRENKIju)

fn <- radni[!(radni$Protein %in% tp$Protein) & radni$Protein %in% frenki$UniProt,]

fn1 <- radni[!(radni$Protein %in% tp$Protein) & 
               radni$Protein %in% frenki_pravo$UniProt,]

## Evaluacija za FRENKI dataset original
precision <- nrow(tp)/(nrow(tp)+nrow(fp))
recall <- nrow(tp)/(nrow(tp)+nrow(fn))
f <- 2*precision*recall/(precision+recall)

## Evaluacija za FRENKI dataset sa izbacenim treningom iz CAFAe
precision1 <- nrow(tp1)/(nrow(tp1)+nrow(fp1))
recall1 <- nrow(tp1)/(nrow(tp1)+nrow(fn1))
f1 <- 2*precision1*recall1/(precision1+recall1)




