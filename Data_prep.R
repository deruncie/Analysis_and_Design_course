
```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(ggplot2)
library(emmeans)
library(lme4)
library(lmerTest)
library(pbkrtest)
library(PLS205)
```

```{r}
data_raw = read.delim('~/Downloads/doi_10_5061_dryad_65d76__v20161202/Experiment1.txt',sep='\t',row.names=1)
```

```{r}
bolting_data = subset(data_raw,Treatment.V %in% c('22ConLDNV','22ConLDV','22VarLDNV','22VarLDV') & Genotype %in% c('Col Ama','Col FRI','flc-3','flc-3 FRI','vin3-4','vin3-4 FRI')) #,'hua2-3','hua2-3 FRI'
bolting_data$Genotype[bolting_data$Genotype == 'Col Ama'] = 'Col'
genotypes = data.frame(Genotype = c('Col','Col FRI','flc-3','flc-3 FRI','vin3-4','vin3-4 FRI'),
                       FRI = c('fri','FRI','fri','FRI','fri','FRI'),
                       mutant = c('WT','WT','flc-3','flc-3','vin3-4','vin3-4'))
bolting_data$FRI = genotypes$FRI[match(bolting_data$Genotype,genotypes$Genotype)]
bolting_data$mutant = genotypes$mutant[match(bolting_data$Genotype,genotypes$Genotype)]
# ggplot(bolting_data,aes(x=Genotype,y=Days.to.Bolt)) +geom_boxplot() + facet_grid(~Vernalization)
bolting_data$Chamber = as.numeric(as.factor(interaction(bolting_data$Chamber.ID,bolting_data$Treatment.V,drop=TRUE)))
cleaned_data = data.frame(Pot = runif(nrow(bolting_data)),bolting_data[,c('Genotype','FRI','mutant','Treatment.V','Chamber','Days.to.Bolt','Chamber.Irradiance')]) #'Temperature','Fluctuation','Vernalization',
cleaned_data = cleaned_data[order(paste(cleaned_data$Vernalization,cleaned_data$Treatment.V,cleaned_data$Chamber),cleaned_data$Pot),]
cleaned_data$Pot = (1:nrow(cleaned_data)-1) %% 24 + 1
write.csv(cleaned_data,file = 'data/Bolting_experiment.csv',row.names=F)
write.csv(subset(cleaned_data,Treatment.V == '22ConLDNV'),file = 'data/Bolting_experiment_1.csv',row.names=F)
write.csv(subset(cleaned_data,Treatment.V == '22ConLDV'),file = 'data/Bolting_experiment_2.csv',row.names=F)

bolting_data_FRI = subset(cleaned_data,Genotype %in% c('Col','Col FRI'))
bolting_data_FRI = bolting_data_FRI[order(bolting_data_FRI$Chamber,bolting_data_FRI$Pot),]
write.csv(bolting_data_FRI,file = 'data/Bolting_experiment_3.csv',row.names=F)

bolting_data = read.csv('data/Bolting_experiment.csv')

mod_NV = lm((Days.to.Bolt)~mutant:FRI,subset(bolting_data,Vernalization == 'NV' & Fluctuation == 'Con'))
means_NV = emmeans(mod_NV,specs = 'FRI',by = 'mutant')
effects = contrast(means_NV,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
regrouped_effects = update(effects,by = 'FRI_effect')
plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)

mod_NV = lm(inverse(Days.to.Bolt)~mutant*FRI,subset(bolting_data,Vernalization == 'NV' & Fluctuation == 'Con'))
means_NV = emmeans(mod_NV,specs = 'FRI',by = 'mutant')
effects = contrast(means_NV,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
regrouped_effects = update(effects,by = 'FRI_effect')
plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)

# mod_NV = lm(log(Days.to.Bolt)~mutant*FRI,subset(bolting_data,Vernalization == 'NV' & Fluctuation == 'Con'))
# means_NV = emmeans(mod_NV,specs = 'FRI',by = 'mutant')
# effects = contrast(means_NV,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
# regrouped_effects = update(effects,by = 'FRI_effect')
# plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)
interactions = contrast(regrouped_effects,'pairwise')
interactions

```
```{r}

mod_V = lm((Days.to.Bolt)~mutant*FRI,subset(bolting_data,Vernalization == 'V' & Fluctuation == 'Con'))
means_V = emmeans(mod_V,specs = 'FRI',by = 'mutant')
effects = contrast(means_V,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
regrouped_effects = update(effects,by = 'FRI_effect')
plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)
mod_V = lm(inverse(Days.to.Bolt)~mutant*FRI,subset(bolting_data,Vernalization == 'V' & Fluctuation == 'Con'))
means_V = emmeans(mod_V,specs = 'FRI',by = 'mutant')
effects = contrast(means_V,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
regrouped_effects = update(effects,by = 'FRI_effect')
plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)
# mod_V = lm(log(Days.to.Bolt)~mutant*FRI,subset(bolting_data,Vernalization == 'V' & Fluctuation == 'Con'))
# means_V = emmeans(mod_V,specs = 'FRI',by = 'mutant')
# effects = contrast(means_V,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
# regrouped_effects = update(effects,by = 'FRI_effect')
# plot(regrouped_effects,horizontal = FALSE,comparisons = TRUE) + geom_vline(xintercept = 0)
interactions = contrast(regrouped_effects,'pairwise')
interactions

```

```{r}
# mod_3way = lmer((Days.to.Bolt)~mutant*FRI*Treatment.V+(1|Chamber.ID:Treatment.V),bolting_data)
# anova(mod_3way,ddf='K')
mod_3way = lmer(inverse(Days.to.Bolt)~mutant*FRI*Treatment.V+(1|Chamber.ID:Treatment.V)+(1|mutant:Chamber.ID:Treatment.V)+(1|FRI:Chamber.ID:Treatment.V),bolting_data)
anova(mod_3way,ddf='K')
# mod_3way = lm((Days.to.Bolt)~mutant*FRI*Treatment.V,bolting_data)
anova(mod_3way,ddf='K')
means = emmeans(mod_3way,specs = 'FRI',by = c('mutant','Treatment.V'))
FRI_effects = contrast(means,'trt.vs.ctrl',ref = 'fri',name = 'FRI_effect')
FRI_effects_regrouped = update(FRI_effects,by = c('FRI_effect','Treatment.V'))
# summary(FRI_effects_regrouped)
mut_effects_on_FRI_effects = contrast(FRI_effects_regrouped,'trt.vs.ctrl',ref = 'FRI - fri WT',by = 'Treatment.V',name = 'mutation_effects')
# summary(mut_effects_on_FRI_effects)
mut_effects_on_FRI_effects_regrouped = update(mut_effects_on_FRI_effects,by = 'mutation_effects')
summary(mut_effects_on_FRI_effects_regrouped)
plot(FRI_effects_regrouped,horizontal = FALSE) + facet_grid(~Treatment.V) + geom_vline(xintercept = 0) + xlab('FRI effect')

plot(mut_effects_on_FRI_effects,horizontal = FALSE) + facet_grid(~Treatment.V) + geom_vline(xintercept = 0) + theme(axis.text.x = element_text(angle=45)) + xlab('mutant effect on FRI effect')
```

```{r}
joint_tests(mod_3way,by = c('Treatment.V'))
joint_tests(mod_3way,by = c('Treatment.V','mutant'))
```


```{r}
data1 = subset(data,Treatment.V %in% c('22VarLDNV') & Genotype %in% c('Col Ama','Col FRI','flc-3','flc-3 FRI','vin3-4','vin3-4 FRI','hua2-3','hua2-3 FRI'))
# data1 = subset(data,Treatment.V %in% c('22VarLDNV') & Genotype %in% c('Col Ama','Col FRI Ama','fca-9','fca-9 FRI','fld-3','fld-3 FRI'))

table(data1$Genotype)
ggplot(data1,aes(x=Genotype,y=1/Days.to.Bolt)) +geom_boxplot()
```

```{r}
data1$FRI = ifelse(grepl('FRI',data1$Genotype),'FRI','fri')
data1$mut = ifelse(grepl('flc-3',data1$Genotype),'flc-3','WT')
data1$mut[grepl('vin3-4',data1$Genotype)] = 'vin3-4'
data1$mut[grepl('hua2-3',data1$Genotype)] = 'hua2-3'
m1 = lm((Days.to.Bolt)~mut*FRI,data1)
means = emmeans(m1,specs = 'FRI',by='mut')
effects = contrast(means,'trt.vs.ctrl',ref='fri',name = 'FRI_effect')
regrouped = update(effects,by = 'FRI_effect')
regrouped
interaction = contrast(regrouped,'trt.vs.ctrl',ref='WT')
interaction
plot(regrouped,comparisons = TRUE,horizontal = FALSE) + geom_vline(xintercept=0)

```
```{r}
means = emmeans(m1,specs = 'mut',by='FRI')
effects = contrast(means,'trt.vs.ctrl',ref='WT',name = 'mut_effect')
regrouped = update(effects,by = 'mut_effect')
regrouped
interaction = contrast(regrouped,'trt.vs.ctrl',ref='fri')
interaction
plot(regrouped,comparisons = TRUE)
```
