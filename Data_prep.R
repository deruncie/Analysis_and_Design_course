library(ggplot2)
library(emmeans)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(PLS205)

data_raw = read.delim('~/Downloads/doi_10_5061_dryad_65d76__v20161202/Experiment1.txt',sep='\t',row.names=1)


bolting_data = subset(data_raw,Treatment.V %in% c('22ConLDNV','22ConLDV','22VarLDNV','22VarLDV') & Genotype %in% c('Col Ama','Col FRI','flc-3','flc-3 FRI','vin3-4','vin3-4 FRI')) #,'hua2-3','hua2-3 FRI'
bolting_data$Genotype[bolting_data$Genotype == 'Col Ama'] = 'Col'
genotypes = data.frame(Genotype = c('Col','Col FRI','flc-3','flc-3 FRI','vin3-4','vin3-4 FRI'),
                       FRI = c('fri','FRI','fri','FRI','fri','FRI'),
                       mutant = c('WT','WT','flc-3','flc-3','vin3-4','vin3-4'))
bolting_data$FRI = genotypes$FRI[match(bolting_data$Genotype,genotypes$Genotype)]
bolting_data$mutant = genotypes$mutant[match(bolting_data$Genotype,genotypes$Genotype)]
# ggplot(bolting_data,aes(x=Genotype,y=Days.to.Bolt)) +geom_boxplot() + facet_grid(~Vernalization)
bolting_data$Chamber = interaction(bolting_data$Chamber.ID,bolting_data$Treatment.V,drop=TRUE)
bolting_data$Chamber = factor(bolting_data$Chamber,levels = sample(levels(bolting_data$Chamber)))
bolting_data$Chamber = as.numeric(bolting_data$Chamber)
cleaned_data = data.frame(Pot = runif(nrow(bolting_data)),bolting_data[,c('Genotype','FRI','mutant','Treatment.V','Chamber','Days.to.Bolt')]) #'Temperature','Fluctuation','Vernalization',
cleaned_data = cleaned_data[order(cleaned_data$Chamber,cleaned_data$Pot),]
cleaned_data$Pot = (1:nrow(cleaned_data)-1) %% 24 + 1
write.csv(cleaned_data,file = 'data/Bolting_experiment_full.csv',row.names=F)
write.csv(subset(cleaned_data,Treatment.V == '22ConLDNV'),file = 'data/Bolting_experiment_1.csv',row.names=F)
write.csv(subset(cleaned_data,Treatment.V == '22ConLDV'),file = 'data/Bolting_experiment_2.csv',row.names=F)
write.csv(subset(cleaned_data,Genotype %in% c('Col','Col FRI')),file = 'data/Bolting_experiment_3.csv',row.names=F)
