## 
##This script reads in a list of files from the directory (the data for each participant in an experiment). 
## The timestamps produced by the participant are converted to inter-event intervals for each dataframe and recurrence quantification analysis 
## is performed on this inter-event-interval sequence. Each participant's timestamps are also randomly shuffled to create
## a control condition and recurrence quantification analysis also applied to this control.
## Finally the list of dataframes is combined into one dataframe 
## and overall trends in the recurrence dynamics for experimental groups are visualised using ggplot plots. 

library(dplyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(crqa)
library(tseriesChaos)
library(entropy)
library(plot3D)
library(lme4)

temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)
## there is now a list of dataframes called myfiles


#Sort the group time stamps in for each trial
ordered_synergy <- lapply(myfiles, function(x) {x$synergytst <- sort(x$synergytst, na.last = TRUE)})
row_nums <- lapply(ordered_synergy, function(x){
  x$rownum <- 1:length(x)
})

# plotting the ordered time stamps to check for misplaced timestamps
lapply(ordered_synergy, function(x) { qplot(1:length(x), x,
                                                  xlab = "clap_number",
                                                  ylab = "point in time") +
    geom_line(size = 0.0001)
  })


#find the inter clap intervals and each trial inter clap interval data as a vector
inter_claps <- lapply(ordered_synergy, function(x){
  #x <- na.omit(x)
  ici <- vector()
  for(i in 1:length(x - 1)){
    ici[i] <- x[i + 1] - x[i]
  }
  x$inter_claps <- ici
  #new_df <- as.data.frame(x, ici)
  #ggplot(new_df, aes(x, ici)) + geom_line()
})

set.seed(1000)

shuffled_claps <- lapply(inter_claps, function(x){
  sample(x)
})

##ARQA function. Set the parameters for all trials here. 
arqa.map <- function(seq) {
  crqa(ts1 = seq, ts2 = seq, delay = 5, embed = 9, rescale = 0, 
       radius = 1, normalize = 0, mindiagline = 2, minvertline = 2,
       tw = 1) 
  # setting tw to 1 removes Line of Identity/Incidence for auto-recurrence analysis
}

## Compute auto recurrence quantification analysis for experimental trial
ARQA_list <- lapply(inter_claps, function(x) {
  arqa.map(na.omit(x))
})

ARQA_list_shuffled <- lapply(shuffled_claps, function(x){
  arqa.map(na.omit(x))
})

## Plot each recurrence plot 
ARQA_Plots <- lapply(ARQA_list, function(x){
  par = list(unit = 3, labelx = "trial time", labely = "trial time", cols = "dark red", pcex = .6, pch = 19, labax = NA, labay = NA, las = 1)
  x.RP <- x$RP
 plot <- plotRP(x.RP, par)
})

ARQA_shuffle_Plots <- lapply(ARQA_list_shuffled, function(x){
  par = list(unit = 3, labelx = "trial time", labely = "trial time", cols = "dark red", pcex = .6, pch = 19, labax = NA, labay = NA, las = 1)
  x.RP <- x$RP
  plot <- plotRP(x.RP, par)
})

## combine the inter clap intervals with the rest of the dataframes
New_Data <- lapply(seq_along(myfiles), function(x){
  cbind(myfiles[[x]], inter_claps[[x]], shuffled_claps[[x]], row_nums[[x]],
        ARQA_list[[x]]$RR,
        ARQA_list[[x]]$DET, 
        ARQA_list[[x]]$LAM,
        ARQA_list[[x]]$ENT,
        ARQA_list[[x]]$maxL,
        ARQA_list_shuffled[[x]]$RR,
        ARQA_list_shuffled[[x]]$DET,
        ARQA_list_shuffled[[x]]$LAM,
        ARQA_list_shuffled[[x]]$ENT,
        ARQA_list_shuffled[[x]]$maxL)
})


## Renaming some variables for convenience
New_Data2 <- bind_rows(New_Data, .id = "column_label")
New_Data2$row_name <- New_Data2$`row_nums[[x]]`
New_Data2$interclapintervals <- New_Data2$`inter_claps[[x]]`
New_Data2$randomclapintervals <- New_Data2$`shuffled_claps[[x]]`

New_Data2$RR <- New_Data2$`ARQA_list[[x]]$RR`
New_Data2$DET <- New_Data2$`ARQA_list[[x]]$DET`
New_Data2$LAM <- New_Data2$`ARQA_list[[x]]$LAM`
New_Data2$ENT <- New_Data2$`ARQA_list[[x]]$ENT`
New_Data2$MaxL <- New_Data2$`ARQA_list[[x]]$maxL`

New_Data2$shuffleRR <- New_Data2$`ARQA_list_shuffled[[x]]$RR`
New_Data2$shuffleDET <- New_Data2$`ARQA_list_shuffled[[x]]$DET`
New_Data2$shuffleLAM <- New_Data2$`ARQA_list_shuffled[[x]]$LAM`
New_Data2$shuffleENT <- New_Data2$`ARQA_list_shuffled[[x]]$ENT`
New_Data2$shuffleMaxL <- New_Data2$`ARQA_list_shuffled[[x]]$maxL`

New_Data2$Patterns <- paste(New_Data2$PatternANum, New_Data2$PatternBNum)

##Plot the inter-event -interval sequence as a facet wrapped line graph.

ggplot(New_Data2, aes(row_name, interclapintervals)) + geom_line(size = 0.0001) + 
  #geom_point(col = "blue", size = 0.0001) + 
  #geom_smooth(method = "lm", size = 0.0001) +
  xlab("clap count") +
  ylab("interclap interval") +
  theme_bw() +
  facet_grid(rows = vars(dyad), cols = vars(trial)) +
  coord_cartesian(ylim = c(0, 3), xlim = c(0, 400))

##Depending on the embedding dimension and delay some Recurrence measures return NAs for some trials,
## na.omit allows us to plot the data anyway, however, this is still an issue we need to deal with.

ggplot(na.omit(New_Data2), aes(condition, na.omit(RR), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
 # geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname) 

ggplot(na.omit(New_Data2), aes(row_name, na.omit(RR), color = condition)) +
  geom_point() +
  theme_bw()

ggplot(na.omit(New_Data2), aes(row_name, na.omit(shuffleRR), color = condition)) +
  geom_point() +
  theme_bw()

ggplot(na.omit(New_Data2), aes(condition, na.omit(shuffleRR), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  # geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname) 

ggplot(na.omit(New_Data2), aes(condition, na.omit(DET), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(shuffleDET), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw()# +
  #facet_wrap(~dyadname)


ggplot(na.omit(New_Data2), aes(condition, na.omit(LAM), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(shuffleLAM), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(ENT), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(shuffleENT), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw()# +
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(MaxL), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() +
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

ggplot(na.omit(New_Data2), aes(condition, na.omit(shuffleMaxL), fill = condition)) +
  geom_jitter(shape = 16, position = position_jitter(0.2)) +
  geom_violin() + 
  #geom_boxplot() +
  theme_bw() #+
  #facet_wrap(~dyadname)

m1 <- glm(RR ~ condition + Patterns,
            data = New_Data2, family = "poisson"(link = log))
summary(m1)

m2 <- glm(shuffleRR ~ condition + Patterns, data = New_Data2,
          family = "poisson"(link = log))
summary(m2)
## Using an optimized generalized mixed effects model to analyse the results.
#m1 <- glmer(ENT ~ condition + Patterns + (1 + condition + Patterns|dyad) +
#              (1 + condition + Patterns|trial),
#            data = New_Data2, family = "poisson"(link = log), glmerControl(optimizer = "bobyqa"))
#summary(m1)
