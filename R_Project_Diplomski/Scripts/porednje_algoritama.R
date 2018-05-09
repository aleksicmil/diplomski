# Poredjenje algoritama

# Ucitava tp i fn iz frenkija i pannzera

f.res <- read.table("Data/frenki_resenja_tpfn.txt", header = T, sep = "\t")
p.res <- read.table("Data/pannzer_resenja_tpfn.txt", header = T, sep = "\t")

# Pravi barplotove koji pokazuju koliko se u tp od frenkija i pannzera nalazi
# poznatih, a koliko predikovanih mitohondrijalnih proteina

barplot(table(f.res[f.res$Guess == "tp", 4]), ylim = c(0, 30), 
        main = "FRENKI - true positives (ukupno 47)")
barplot(table(p.res[p.res$Guess == "tp", 4]), ylim = c(0, 200), 
        main = "PANNZER - true positives (ukupno 230)")

# Testira da li je razika u broju poznatih i predvidjenih znacajna

f.chi <- table(f.res[f.res$Guess == "tp", 4])
# Known mitochondrial Predicted mitochondrial 
# 26                      21 

p.chi <- table(p.res[p.res$Guess == "tp", 4])
# Known mitochondrial Predicted mitochondrial 
# 181                      49

chisq.test(rbind(f.chi, p.chi))
# Pearson's Chi-squared test with Yates' continuity correction
# 
# data:  rbind(f.chi, p.chi)
# X-squared = 10.089, df = 1, p-value = 0.001492
