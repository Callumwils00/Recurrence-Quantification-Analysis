library(tidyverse)
library(purrr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(crqa)
library(tseriesChaos)
library(entropy)
library(plot3D)
library(lme4)
library(cowplot)
library(gridExtra)
library(reshape2)
library(ggpubr)

Top_folder <- "Z:\\Synergy files\\Purrr solution"
Prop_data <- "Z:\\Synergy files\\Propdata"
Recurrence_data <- "Z:\\Synergy files"
output_folder <- "\\output"
synergyfiles <- "\\syn_out\\"
Autofiles <- "\\ARQA_out\\"
CRQAfiles <- "\\CRQA_out\\"
SubDatafiles <- "\\SubData_out\\"

dir.create(file.path(Prop_data), showWarnings = FALSE)
dir.create(file.path(Recurrence_data), showWarnings = FALSE)
dir.create(file.path(Top_folder, output_folder), showWarnings = FALSE)
dir.create(file.path(Top_folder, synergyfiles), showWarnings = FALSE)
dir.create(file.path(Top_folder, Autofiles), showWarnings = FALSE)
dir.create(file.path(Top_folder, CRQAfiles), showWarnings = FALSE)
dir.create(file.path(Top_folder, SubDatafiles), showWarnings = FALSE)

Datasynch <- read.csv(file = file.path(Prop_data, "claprevisedSimpleComplex.csv"))
Datasynch$dyadtrial <- paste0(Datasynch$dyad, "_", Datasynch$trial)
Datasynchlist <- Datasynch %>% split(.$dyadtrial)

proportionate <- function(x){
  ##creating a clap count variable and a time variable from the first clap in the trial to the last clap.
  if(max(na.omit(x$PartnerA)) >= max(na.omit(x$PartnerB))){
    clapcount <- seq(1:length(x$PartnerA)) 
    x$clapcount <- clapcount
    max <- max(na.omit(x$PartnerA))
  } else {
    clapcount <- seq(1: length(x$PartnerB))
    x$clapcount <- clapcount
    max <- max(na.omit(x$PartnerB))
  }
  
  if(min(na.omit(x$PartnerA)) <= min(na.omit(x$PartnerB))){
    min <- min(na.omit(x$PartnerA))
  } else {
    min <- min(na.omit(x$PartnerB))
  }
  ##the time variable
  time <- seq(min, max, length.out = length(clapcount))
  ##a max clap variable, as now the absolute number of claps is a potential random effect whe using bin number
  max_claps <- max(clapcount)
  ## creating 50 uniformly spaced bins in each trial
  proportion_bins = ntile(clapcount, 10)
  
  x$proportion_bins <- as.factor(proportion_bins)
  ##create synch to simple and synch to complex for each bin
  x <- x %>% group_by(proportion_bins) %>% mutate(propsynchsimple = mean(SyncSimple == 1)*100,
                                                  propsynchcomplex = mean(SyncComplex == 1)*100)
  
  x <- as.data.frame(x)
  x <- cbind(x, time)
  x <- cbind(x, max_claps)
  x <- x[!duplicated(x$proportion_bins), ]
}

## bind all the trials together into a new full dataframe 
newDatasynch <- map_dfr(Datasynchlist, proportionate)

##Plotting to viausalize trends

colors <- c("proportion sync to simple" = "blue", "proportion sync to complex" = "red")

newDatasynch$propsynchsimple <- as.numeric(newDatasynch$propsynchsimple)

ggplot(newDatasynch, aes(proportion_bins, propsynchsimple)) + 
  geom_point(stat='summary', fun = mean, aes(color = "proportion sync to simple")
  ) +
  geom_line(stat='summary', fun = mean, aes(color = "proportion sync to simple", group = 1)
  ) +
  geom_point(aes(proportion_bins, propsynchcomplex, color = "proportion sync to complex"),
             stat = 'summary', fun= mean) +
  geom_line(aes(proportion_bins, propsynchcomplex, color = "proportion sync to complex", group = 1),
            stat = 'summary', fun= mean) +
  geom_smooth(aes(group = 1),method = "gam", formula = y ~ poly(x, 2)) +
  geom_smooth(aes(proportion_bins, propsynchcomplex, color = "proportion sync to complex", group = 1), 
              method = "gam", formula = y ~ poly(x, 2)) +
  ylab("proportion synch") +
  facet_grid(~condition)


##Change the data types
newDatasynch$proportion_bins <- as.numeric(newDatasynch$proportion_bins)
newDatasynch$propsynchsimple <- as.numeric(newDatasynch$propsynchsimple)
newDatasynch$propsynchcomplex <- as.numeric(newDatasynch$propsynchcomplex)
newDatasynch$condition <- as.factor(newDatasynch$condition)

##Growth Curve analysis models
m <- lmer(propsynchsimple ~ poly(proportion_bins, 2)*condition +
            (1|dyad) + (0 + poly(proportion_bins, 2)|dyad)
          + (1|max_claps) + (0 + poly(proportion_bins,2)|max_claps), data = newDatasynch)
summary(m)

m2 <- lmer(propsynchcomplex ~poly(proportion_bins, 2)*condition +
             (1| dyad) + (0 + poly(proportion_bins, 2)|dyad) +
             (1|max_claps) + (0 + poly(proportion_bins,2)|max_claps),data = newDatasynch)
summary(m2)

#Choose how to import the data
#1: there are many files in the directory
Data <- read.csv(file = file.path(Top_folder, "Anon Full Recurrence dataset2.csv"))
Data$dyadtrial <- paste(Data$dyadname, Data$trial)
myfiles <- Data %>% split(.$dyadtrial)
##there is now a list of dataframes called myfiles

#Set the arqa params for all trials here
arqa.map <- function(seq, de, em) {
  crqa(ts1 = seq, ts2 = seq, delay = de, embed = em, rescale = 0, 
       radius = .8, normalize = 0, mindiagline = 2, minvertline = 2,
       tw = 1) 
  # setting tw to 1 removes Line of Identity/Incidence for auto-recurrence analysis
}

de <- c(1,2)
em <- c(8,10)

#Creating two arqa.map functions allow you to set different parameters for the 
#synergy and auto recurrence calculations 
arqa.map2 <- function(seq, Ade, Aem) {
  crqa(ts1 = seq, ts2 = seq, delay = Ade, embed = Aem, rescale = 0, 
       radius = .8, normalize = 0, mindiagline = 2, minvertline = 2,
       tw = 1) 
  # setting tw to 1 removes Line of Identity/Incidence for auto-recurrence analysis
}

Ade <- c(1, 2)
Aem <- c(8, 10)

#Set the crqa params for all trials here
crqa.map <- function(seq1, seq2, Cde, Cem) {
  crqa(ts1 = seq1, ts2 = seq2, delay = 1, embed = 2, rescale = 1, 
       radius = 1.2, normalize = 1, mindiagline = 2, minvertline = 2,
       tw = 0)
}

Cde <- c(1,1)
Cem <- c(2,2)

##A function that plots the RP for any rqa results passed to it
param = list(unit = 3, labelx = "", labely = "", cols = "dark green",
             pcex =.01, pch = 19, labax = NA, labay = NA, las = 1)


RQAplot <- function(x, y){
  rqa_df = data.frame(points = y$RP@i,
                      loc = seq_along(y$RP@i))
  
  ggplot(rqa_df,aes(x=points,
                    y=loc)) +
    geom_point(color="red",size=.01) +
    ggtitle(toString(unique(x$dyadtrial))) + coord_fixed() +
    theme_classic() +
    theme(legend.position="none", axis.text.x = element_blank(), axis.text.y = element_blank(),
          aspect.ratio = 1)
}

analysis <- function(x) {
  newx <- x %>% select("rownumber", "PartnerA", "PartnerB")
  subdata <- melt(newx, id.vars = c("rownumber"), variable.names = "jointclapping")
  joint_clapping <- subdata$variable
  synergy <- as.numeric(subdata$value)
  trial <- toString(x$dyadtrial)
  sorted_synergy <- sort(synergy, na.last = TRUE)
  
  #Creating each inter clap interval series before analysis
  ici <- vector()
  for(i in 1:length(sorted_synergy - 1)){
    ici[i] <- sorted_synergy[i + 1] - sorted_synergy[i]}
  
  clapA <- x$PartnerA
  clapB <- x$PartnerB
  
  iciA <- vector()
  for(i in 1:length(clapA - 1)){
    iciA[i] <- clapA[i + 1] - clapA[i]}
  
  iciB <- vector()
  for(i in 1:length(clapB - 1)){
    iciB[i] <- clapB[i + 1] - clapB[i]}
  
  iciA2 <- vector()
  for(i in 1:length(clapA - 1)){
    iciA2[i] <- clapA[i + 1] - clapA[i]}
  
  iciB2 <- vector()
  for(i in 1:length(clapB - 1)){
    iciB2[i] <- clapB[i + 1] - clapB[i]
  }
  
  iciA2[is.na(iciA2)] <- 50
  iciB2[is.na(iciB2)] <- -50
  
  for(i in de){
    for(p in em){
      synName <- paste0("syn.arqa", "_", toString(i), "_", toString(p), "_")
      synNameRes <- arqa.map(ici, i, p)
      RR <- synNameRes$RR
      DET <- synNameRes$DET
      
      new_df <- cbind(sorted_synergy,
                      ici, 
                      joint_clapping,
                      RR, DET, toString(i),
                      toString(p),
                      paste0(toString(i), "_", toString(p)))
      
      colnames(new_df) <- c("synergy","ici",
                            "joint_clapping","syn_arqa.RR", "syn_arqa.DET",
                            "delay", "embed", "del_emb"
      )
      
      new_df <- cbind(x, new_df)
      write.csv(new_df, file= paste0(Top_folder, synergyfiles,
                                     paste0(synName,toString(unique(x$dyadtrial)),".csv"))) ## Sending each analysis to a separate folder,
    } 
  }
  
  for(i in Ade){
    for(p in Aem){
      partners <- list(iciA, iciB)
      partnernames <- c("iciA", "iciB")
      for(q in partners){
        for(w in partnernames){
          ARQAName <- paste0(w, "_", toString(i), "_", toString(p), "_")
          ARQARes <- arqa.map2(na.omit(q), i, p)
          RR <- ARQARes$RR
          DET <- ARQARes$DET
          
          new_df <- cbind(sorted_synergy,
                          iciA, 
                          joint_clapping,
                          RR, DET, toString(i),
                          toString(p),
                          paste0(toString(i), "_", toString(p)))
          
          colnames(new_df) <- c("synergy","iciA",
                                "joint_clapping","ARQA.RR", "ARQA.DET",
                                "delay", "embed", "del_emb"
          )
          
          new_df <- cbind(x, new_df)
          
          write.csv(new_df, file= paste0(Top_folder, Autofiles,
                                         paste0(ARQAName,toString(unique(x$dyadtrial)),".csv")))
        }
      }
    }
  }
  
  for(i in Cde){
    for(p in Cem){
      synchName <- paste0("synch.crqa", "_", toString(i), "_", toString(p), "_")
      synchNameRes <- crqa.map(iciA2, iciB2, i, p)
      RR <- synchNameRes$RR
      DET <- synchNameRes$DET
      
      new_df <- cbind(sorted_synergy,
                      ici, 
                      joint_clapping,
                      RR, DET, toString(i),
                      toString(p),
                      paste0(toString(i), "_", toString(p)))
      
      colnames(new_df) <- c("synergy","ici",
                            "joint_clapping","syn_arqa.RR", "syn_arqa.DET",
                            "delay", "embed", "del_emb"
      )
      
      new_df <- cbind(x, new_df)
      write.csv(new_df, file= paste0(Top_folder, CRQAfiles,
                                     paste0(synchName,toString(unique(x$dyadtrial)),".csv"))) ## Sending each analysis to a separate folder,
    } 
  }
  
  ## Below is just creating some plots to visualize the data for each dyad_trial
  syn.arqa <- arqa.map(ici, 1, 10)
  subdata$ici <- ici
  P1 <- RQAplot(x, syn.arqa)
  subdata$clap_num1 <- 1:length(ici)
  PL1 <- ggplot(subdata, aes(clap_num1, ici)) + geom_line() +
    xlab("clap number") +
    ylab("Inter claps") +
    theme_bw()+
    ggtitle("Synergy Series")
  
  partA.arqa <- arqa.map2(na.omit(iciA), 1, 10)
  P2 <- RQAplot(x, partA.arqa)
  subdata$clap_numA <- 1:length(iciA)
  subdata$iciA <- iciA
  PL2 <-  ggplot(subdata, aes(clap_numA, iciA)) + geom_line() +
    xlab("clap number") +
    ylab("Inter claps") +
    theme_bw()+
    ggtitle("Partner A Series")
  
  partB.arqa <- arqa.map2(na.omit(iciB), 1, 10)
  P3 <- RQAplot(x, partB.arqa)
  subdata$clap_numB <- 1:length(iciB)
  subdata$iciB <- iciB
  P3L <-  ggplot(subdata, aes(clap_numB, iciB)) + geom_line() +
    xlab("clap number") +
    ylab("Inter claps") +
    theme_bw()+
    ggtitle("Partner B Series")
  
  iciA2 <- vector()
  for(i in 1:length(clapA - 1)){
    iciA2[i] <- clapA[i + 1] - clapA[i]}
  
  iciB2 <- vector()
  for(i in 1:length(clapB - 1)){
    iciB2[i] <- clapB[i + 1] - clapB[i]
  }
  
  iciA2[is.na(iciA2)] <- 50
  iciB2[is.na(iciB2)] <- -50
  
  synch.crqa <- crqa.map(iciA2, iciB2)
  P4 <- RQAplot(x, synch.crqa)
  PL4 <- ggplot(subdata) + geom_line(aes(x = clap_numA, y = iciA), color = "blue") +
    geom_line(aes(x = clap_numA, y = iciB), color = "red") +
    xlab("clap number") +
    ylab("Inter claps") +
    theme_bw()+
    ggtitle("Synchrony CRQA")
  
  Panel <- ggdraw() +
    draw_plot(P1,  x = 0, y = .8, width = .2, height = .2) +
    draw_plot(PL1, x = .2, y = .8, width = .2, height = .2) +
    draw_plot(P2, x = 0, y = .6, width = .2, height = .2) +
    draw_plot(PL2, x = .2, y = .6, width = .2, height = .2) +
    draw_plot(P3, x = 0, y = .4, width =.2, height = .2) +
    draw_plot(P3L, x = .2, y = .4, width =.2, height = .2) +
    draw_plot(P4, x =0, y = .2, width = .2, height = .2) +
    draw_plot(PL4, x =.2, y =.2, width =.2, height = .2)
  print(Panel, width = 1920/72, height = 1080/72, dpi = 72)
  ggsave(file = paste0("Panel", "_", toString(unique(x$dyadtrial)), ".png"), plot = Panel,
         path = paste0(Top_folder, output_folder),
         device = "png",
         width = 1920/72, height = 1080/72, dpi = 72)
}

new_data <- purrr::map(myfiles, analysis)






##Reading in all the synergy analysis files
tempsyn <- list.files(path = paste0(Top_folder, synergyfiles), pattern="*.csv")
mysynergyfiles <- lapply(tempsyn, function(x){
  files <- paste0(Top_folder, synergyfiles, x)
  read.csv(files)
})

##Create on df
synergydf <-bind_rows(mysynergyfiles, .id = "column_label")

## Look at the correlation between Recurrence measures
ggplot(synergydf, aes(syn_arqa.RR, syn_arqa.DET, color = condition)) + geom_point() + geom_smooth(method = "glm") +
  xlab("Recurrence Rate") +
  ylab("Determinism") +
  ggtitle("Visualisation with Conditions Synergy") +
  theme_bw()

ggplot(synergydf, aes(syn_arqa.RR, syn_arqa.DET, color = del_emb)) + geom_point() + geom_smooth(method = "glm") +
  xlab("Recurrence Rate") +
  ylab("Determinism") +
  ggtitle("Visualisation with Delay and Embedding Dim Synergy") +
  theme_bw()

##Split by delay and embedding dim combination. To create a list of del_emb defined dfs.
new_synergydf <- synergydf %>% split(.$del_emb)

## We can define the regression model in a function and then apply it to the data under each combination of embed and delay.

##Same for auto recurrence analysis
tempAuto <- list.files(path = paste0(Top_folder, Autofiles), pattern="*.csv")
myautofiles <- lapply(tempAuto, function(x){
  files <- paste0(Top_folder, Autofiles, x)
  read.csv(files)
})

Autofilesdf <-bind_rows(myautofiles, .id = "column_label")
## Look at the correlation between Recurrence measures
ggplot(Autofilesdf, aes(ARQA.RR, ARQA.DET, color = condition)) + geom_point() + geom_smooth(method = "glm") +
  xlab("Recurrence Rate") +
  ylab("Determinism") +
  ggtitle("Visualisation with Conditions Auto Recurrence") +
  theme_bw()


ggplot(Autofilesdf, aes(ARQA.RR, ARQA.DET, color = del_emb)) + geom_point() + geom_smooth(method = "glm") +
  xlab("Recurrence Rate") +
  ylab("Determinism") +
  ggtitle("Visualisation with Delay and Embedding Dim Auto Recurrence") +
  theme_bw()
