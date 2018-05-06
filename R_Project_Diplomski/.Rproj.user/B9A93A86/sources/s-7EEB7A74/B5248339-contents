
###_________________________________________________________
### Mapiranje za evaluacije
###_________________________________________________________

mapa <- read.table("Data/CAFA_mapping.txt", head <- F, sep = "\t", 
                   col.names = c("UniProt", "Uniprot1", "UniProt2", "CAFA"))

frenki <- read.table("Data/FRENKI_CC_predikcije_propagirano.txt", header = F, 
                     sep= "\t", col.names = c("CAFA", "Predikcija", "Skor"))

# Frenki je anotiran po CAFA anotacijama. fajl mapa poveyuje CAFA anptacije i
# uniprot koje meni trebaju. Treba frenki fajlu da dodam uniprot anotacije.

frenki_score <- merge(frenki, mapa, by = "CAFA", all.x = T)
write.table(frenki_score, "Data/FRENKI_predikcije_mapirano.txt", sep = "\t", 
            row.names = F, col.names = T, quote = F)
