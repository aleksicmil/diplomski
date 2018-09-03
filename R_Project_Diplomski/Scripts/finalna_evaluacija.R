###
### FINALNA EVALUACIJA ALATA 
###
###

# ________________________________________________________________________
# UCITAVANJE DATASETOVA
# ________________________________________________________________________

# Referentni dataset mitohondrijalnih proteina. Napravljen iz liste IMPI gena 
# prevedene u odgovarajuce proteine, kojoj su oduzeti svi proteini koji su
# trenutno u GO anotirani sa GO:0005739 ili se nalaze u CAFA trening setu.

ref <- read.table(ref, "Data/FINALNO/referentni_mito_finalno.txt", 
                  row.names = F, col.names = T, quote = F, sep = "\t")

# Finalni dataset koji sadrzi FRENKI predikcije. Predikcije su propagirane do
# roditeljskih termina, i ekstrahovane su samo predickije za mitohondrijalnu
# lokalizaciju. Izbaceni su duplikati, NA, svi proteini iz CAFA trening seta, 
# i svi proteini u GO trenutno anotirani sa GO:0005739.

frenki <- read.table("Data/FINALNO/frenki_pred_finalno.txt", header = T, 
                     sep = "\t")

# Finalni dataset koji sadrzi PANNZER predikcije. Predikcije su propagirane do
# roditeljskih termina, i ekstrahovane su samo predickije za mitohondrijalnu
# lokalizaciju. Izbaceni su duplikati, NA, svi proteini iz CAFA trening seta, 
# i svi proteini u GO trenutno anotirani sa GO:0005739.

pan <- read.table("Data/FINALNO/pan_pred_fin.txt", header = T, sep = "\t")

# ________________________________________________________________________
# EVALUACIJA
# ________________________________________________________________________


## FRENKI 

#  True positives (99) - tacne predikcije
f.tp <- as.data.frame(intersect(frenki$UniProt, ref$UniProt)) 

#  False positives (2648) - pogresne predikcije
f.fp <- as.data.frame(setdiff(frenki$UniProt, ref$UniProt)) 

# False negatives (554) - proteini koji nisu predvidjeni kao mito, a u 
# referentnom datasetu su
f.fn <- as.data.frame(setdiff(ref$UniProt, frenki$UniProt)) 

# Preciznost (0.036) - prec = tp/(tp+fp)
f.prec <- 99/(99+2648) #0.036

# Odziv (0.151) - recall = tp/(tp+fn)
f.rec <- 99/(99+554) #0.151

# F mera (0.058) - f=2*prec*rec/(prec+rec)
f.f <- 2*f.prec*f.rec/(f.prec+f.rec) #0.058


## PANNZER


#  True positives (343) - tacne predikcije
p.tp <- as.data.frame(intersect(pan$UniProt, ref$UniProt)) 

#  False positives (608) - pogresne predikcije
p.fp <- as.data.frame(setdiff(pan$UniProt, ref$UniProt)) 

# False negatives (310) - proteini koji nisu predvidjeni kao mito, a u 
# referentnom datasetu su
p.fn <- as.data.frame(setdiff(ref$UniProt, pan$UniProt)) 

# Preciznost (0.360) - prec = tp/(tp+fp)
p.prec <- 343/(343+608) 

# Odziv (0.525) - recall = tp/(tp+fn)
p.rec <- 343/(343+310) 

# F mera (0.427) - f=2*prec*rec/(prec+rec)
p.f <- 2*p.prec*p.rec/(p.prec+p.rec) #0.427

#_________________________________________________________________________
# PREKLAPANJA
#_________________________________________________________________________


# Preklapanje izmedju true positive-a FRENKI-ja i PANNZER-a
tp.oba <- as.data.frame(intersect(p.tp[, 1], f.tp[ , 1]))

# Preklapanje izmedju false positive-a FRENKI-ja i PANNZER-a
fp.oba <- as.data.frame(intersect(p.fp[, 1], f.fp[ , 1]))

# Proteini koje nisu predikovali ni PANNZER ni FRENKI a u referentnom 
# datasetu su
no.guess <- as.data.frame(setdiff(ref$UniProt, pan$UniProt))
names(no.guess) <- "UniProt"
no.guess <- as.data.frame(setdiff(no.guess$UniProt, frenki$UniProt))
names(no.guess) <- "UniProt"
no.guess <- merge(no.guess, impi, by = "UniProt", all.x = T )

write.table(no.guess$Gene, "Data/FINALNO/no_guess_gene.txt", sep = "\t", 
            col.names = F, row.names = F, quote = F)

