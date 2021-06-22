library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)
# install.packages("ggpubr")
setwd("~/PycharmProjects/StageORFphase")
tab_kmer_all <- data.frame()
wt_phasing <- data.frame()
wt_nc <- data.frame()
SRR_tab = read.csv("TrueAdapt.csv", header=T)
# wt_phasing <- read.table("SRR1520311_phase_median_mean.tab", sep = "\t", header = T)
# wt_nc <- read.table("./kmer_28/SRR1520311_kmer_28_nc_reads.tab", header= T, sep = "\t")

for (i in c("SRR1520311","SRR1520312", "SRR1520313", "SRR1520314", "SRR1520315", "SRR1520316", "SRR1520317", "SRR1520318")){
  tab_kmer <- read.table(paste(i,"_phase_median_mean.tab",sep=""),sep = "\t", header=T)
  size_nc <- tab_kmer[which.max(tab_kmer$Median_p0),"Size"]
  tab_phasing <- read.table(paste(i,"_phase_median_mean.tab",sep = ""), sep = "\t", header = T)
  tab_phasing["Type"] <-  "WT"
  tab_nc <- read.table(paste("./kmer_",size_nc,"/",i,"_kmer_",size_nc,"_nc_reads.tab",sep = ""), header= T, sep = "\t")
  tab_nc["SRR"] <- i
  tab_nc["Type"] <- "WT"
  wt_nc <- rbind(wt_nc,tab_nc)
  wt_phasing <- rbind(wt_phasing,tab_phasing)
  tab_kmer_all <- rbind(tab_kmer_all,tab_kmer)
}

stress_phasing <- data.frame()
stress_nc <- data.frame()
for (i in c("SRR1520319","SRR1520320","SRR1520321","SRR1520322","SRR1520323","SRR1520324","SRR1520325","SRR1520326","SRR1520327","SRR1520328")){
  tab_kmer <- read.table(paste(i,"_phase_median_mean.tab",sep=""),sep = "\t", header=T)
  size_nc <- tab_kmer[which.max(tab_kmer$Median_p0),"Size"]
  tab_phasing <- read.table(paste(i,"_phase_median_mean.tab",sep = ""), sep = "\t", header = T)
  tab_phasing["Type"] = "Stress"
  tab_nc <- read.table(paste("./kmer_",size_nc,"/",i,"_kmer_",size_nc,"_nc_reads.tab",sep = ""), header= T, sep = "\t")
  tab_nc["SRR"] <- i
  tab_nc["Type"] <- "Stress"
  stress_nc <- rbind(stress_nc,tab_nc)
  stress_phasing <- rbind(stress_phasing,tab_phasing)
  tab_kmer_all <- rbind(tab_kmer_all,tab_kmer)
}


wt_stress_phasing = rbind(wt_phasing,stress_phasing)
wt_stress_phasing = pivot_longer(wt_stress_phasing,cols=c(starts_with("Mean"),starts_with("Median")), names_to = c("Methode","Phase"), names_sep = c("_"), values_to = "Value")
wt_stress_nc = rbind(wt_nc,stress_nc)

ggplot(wt_stress_phasing[wt_stress_phasing["Methode"]=="Mean",], aes(x=Size, y=Value, fill=Riboseq))+
  geom_bar(position="dodge", stat = "identity")+
  facet_grid(vars(Phase), vars(Type))+
  ggtitle("Mean per size per phase per type")
ggplot(wt_stress_phasing[wt_stress_phasing["Methode"]=="Median",], aes(x=Size, y=Value, fill=Riboseq))+
  geom_bar(position="dodge", stat = "identity")+
  facet_grid(vars(Phase), vars(Type))+
  ggtitle("Median per size per phase per type")

good_nc_number = function(number, tab=wt_stress_nc){
  good_nc = tab[tab["Number.reads"]>=number,]
  # wt_good_nc_10 <-pivot_longer(wt_good_nc_10,cols = contains(".p"), names_to = c("Perc.","Number"), names_sep = ".", values_to = "Perc.")
  good_nc <-pivot_longer(good_nc,cols = c("Perc..p0","Perc..p1","Perc..p2"), names_to = "Phase", names_prefix = ("Perc.."), values_to = "Perc.")
  good_nc <-pivot_longer(good_nc,cols = c("Number.p0","Number.p1","Number.p2"), names_to = "Phase1", names_prefix = ("Number."), values_to = "Number.")
  good_nc <- good_nc[which(good_nc$Phase == good_nc$Phase1),]
  good_nc <-select(good_nc,-Phase1) 
  return(good_nc)
}
good_nc_number(100)
# wt_good_nc_10 = wt_nc[wt_nc["Number.reads"]>=10,]
# # wt_good_nc_10 <-pivot_longer(wt_good_nc_10,cols = contains(".p"), names_to = c("Perc.","Number"), names_sep = ".", values_to = "Perc.")
# wt_good_nc_10 <-pivot_longer(wt_good_nc_10,cols = c("Perc..p0","Perc..p1","Perc..p2"), names_to = "Phase", names_prefix = ("Perc.."), values_to = "Perc.")
# wt_good_nc_10 <-pivot_longer(wt_good_nc_10,cols = c("Number.p0","Number.p1","Number.p2"), names_to = "Phase1", names_prefix = ("Number."), values_to = "Number.")
# wt_good_nc_10 <- wt_good_nc_10[which(wt_good_nc_10$Phase == wt_good_nc_10$Phase1),]
# wt_good_nc_10 <-select(wt_good_nc_10,-Phase1) 
str(good_nc_number(10)[wt_good_nc_number(10)$Phase=="p0","Perc."])
nc <- ggplot(good_nc_number(0),aes(x=Perc.,fill = Phase))+
  geom_density()+
  facet_grid(vars(Type))
phasing <- ggplot(good_nc_number(100),aes(x=Perc.,fill = Phase))+
  geom_density()
ggplot(wt_good_nc_number(50),aes(y=Perc., fill = Phase))+
  geom_boxplot()
wt_good_nc_10_90 = wt_nc[wt_nc["Number.reads"]>=10 & wt_nc["Perc..p0"]>=90,]
wt_good_nc_10_80 = wt_nc[wt_nc["Number.reads"]>=10 & wt_nc["Perc..p0"]>=80,]
wt_good_nc_30_80 = wt_nc[wt_nc["Number.reads"]>=30 & wt_nc["Perc..p0"]>=80,]
wt_good_nc <-gather(wt_good_nc,key = "Phase", value = "Percentage", Perc..p0,Perc..p1,Perc..p2)
wt_nc <-gather(wt_nc,key = "Phase", value = "Percentage", Perc..p0,Perc..p1,Perc..p2)
ggplot(wt_nc, aes(x=Number.reads))+
  geom_histogram(breaks =seq(10,50,by=2))
ggplot(wt_nc, aes(x=Number.reads))+
  geom_histogram(breaks =seq(1,50,by=2))
table(wt_nc["Number.reads"])
ggplot(wt_nc[wt_nc["Number.reads"]>=10,],aes(x=Phase, y=Percentage))+
  geom_boxplot()
barplot(table(wt_nc["Number.reads"]))
ggplot(wt_good_nc,aes(x=Phase,y=Percentage))+
  geom_boxplot()
wt_nc_agg = aggregate(ID ~ Number.reads, data = wt_coding, FUN = length)


polyApO <- polyA[polyA$Phase == 0,]
polyAp1 <- polyA[polyA$Phase == 1,]
polyAp2 <- polyA[polyA$Phase == 2,]
polyA$Phase <- as.factor(polyA$Phase)
pmed <- ggplot(polyA, aes(x=Phase,y=Median))+
  geom_boxplot()
pmed
pmean <-ggplot(polyA, aes(x=Phase,y=Mean))+
  geom_boxplot()
pmean

polyA_all <- read.csv("best_phasemean_p0_reads_polyA.csv", header = TRUE)
pperc <- ggplot(polyA_all, aes(y=Percentage.of.p0))+
  geom_boxplot()
pperc

set1 <- read.table()
