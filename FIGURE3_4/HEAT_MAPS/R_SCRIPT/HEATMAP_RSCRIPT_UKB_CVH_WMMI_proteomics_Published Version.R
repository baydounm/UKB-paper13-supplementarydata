#######################################################################
#######  Making heatmaps for UKB project #13 CVH and WM microstructure integrity mediated by proteomics
#######  Mission: Plot the betas and mediation effects
#######  Programmer: Yi-Han Hu
#######  Date: Mar. 14 2024
#######  revised based on May's suggestions: Mar. 18 2024
#######  updated after updating the LE8 score: Aug. 3 2024
#######  updated for response - all %mediated become solid if TE pvalue <0.05: Oct. 18 2024
#######################################################################

op <- options(nwarnings = 10000)
# --------------------------------------
# Specify working directory where the script and data files are
# --------------------------------------
WorkingDirectory = "file route"

# --------------------------------------
# Set working directory
# --------------------------------------
setwd(WorkingDirectory)

# --------------------------------------
# simple function
# --------------------------------------
floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)

# --------------------------------------
# Turn off scientific notation
# --------------------------------------
options(scipen=999)

# --------------------------------------
# Install/load the packages
# --------------------------------------
library(haven)
library(tidyr) 
library(dplyr)
library(stringr)
library(ggplot2)
library(reshape)
library(purrr)
library(data.table)
library(sjmisc)
library(RColorBrewer)
library(ggnewscale)
library(scico)
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------- Part 1 Data preprocess -----------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# --------------------------------------
# Load data:
# --------------------------------------

## 1.	No pure mediation: non-significant PIE
# NO_MEDIATION_FA <- read_dta("Data/NO_MEDIATION_GROUPA_FAwide.dta")
NO_MEDIATION_FA <- read_dta("UKB_paper13_heatmap_updating/NO_MEDIATION_GROUPA_FAwide.dta")

dim(NO_MEDIATION_FA)
head(NO_MEDIATION_FA)
colnames(NO_MEDIATION_FA)

NO_MEDIATION_FA <- NO_MEDIATION_FA %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

NO_MEDIATION_FA.long <- to_long(NO_MEDIATION_FA, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE'))

## 2.	Inconsistent mediation: PIE is significant but CDE>TE because PIE goes in the opposite direction to TE.
# FA
# INCONSISTENT_MEDIATION_FA <- read_dta("Data/INCONSISTENT_MEDIATION_GROUPB_FAwide.dta")
INCONSISTENT_MEDIATION_FA <- read_dta("UKB_paper13_heatmap_updating/INCONSISTENT_MEDIATION_GROUPB_FAwide.dta")
dim(INCONSISTENT_MEDIATION_FA)
head(INCONSISTENT_MEDIATION_FA)
colnames(INCONSISTENT_MEDIATION_FA)


INCONSISTENT_MEDIATION_FA <- INCONSISTENT_MEDIATION_FA %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

INCONSISTENT_MEDIATION_FA.long <- to_long(INCONSISTENT_MEDIATION_FA, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE'))

# OD
# INCONSISTENT_MEDIATION_OD <- read_dta("Data/INCONSISTENT_MEDIATION_GROUPAB_ODwide.dta")
INCONSISTENT_MEDIATION_OD <- read_dta("UKB_paper13_heatmap_updating/INCONSISTENT_MEDIATION_GROUPAB_ODwide.dta")
dim(INCONSISTENT_MEDIATION_OD)
head(INCONSISTENT_MEDIATION_OD)
colnames(INCONSISTENT_MEDIATION_OD)


INCONSISTENT_MEDIATION_OD <- INCONSISTENT_MEDIATION_OD %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

INCONSISTENT_MEDIATION_OD.long <- to_long(INCONSISTENT_MEDIATION_OD, keys = 'term',
                                          values = c('estimate','p'), 
                                          c('TE','CDE','INTREF', 'INTMED', 'PIE'),
                                          c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE'))

## 3.	Consistent mediation: PIE is significant and CDE<TE because PIE goes in the same direction as TE
# CONSISTENT_MEDIATION_FA <- read_dta("Data/CONSISTENT_MEDIATION_GROUPC_FAwide.dta")
CONSISTENT_MEDIATION_FA <- read_dta("UKB_paper13_heatmap_updating/CONSISTENT_MEDIATION_GROUPC_FAwide.dta")
dim(CONSISTENT_MEDIATION_FA)
head(CONSISTENT_MEDIATION_FA)
colnames(CONSISTENT_MEDIATION_FA)

CONSISTENT_MEDIATION_FA <- CONSISTENT_MEDIATION_FA %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie, betap_intref, betap_pie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie, percent_interaction = betap_intref, percent_mediated = betap_pie) %>% 
  mutate(p_mediated = ifelse(!is.na(p_TE) & p_TE < 0.05, 0.01, 0.1),
         p_interaction = ifelse(!is.na(p_TE) & p_TE < 0.05 & !is.na(p_INT) & p_INT < 0.05, 0.01, 0.1)) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

CONSISTENT_MEDIATION_FA.long <- to_long(CONSISTENT_MEDIATION_FA, keys = 'term',
                                  values = c('estimate','p'), 
                                  c('TE','CDE','INTREF', 'INTMED', 'PIE', 'percent_interaction', 'percent_mediated'),
                                  c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE', 'p_interaction', 'p_mediated'))


# CONSISTENT_MEDIATION_OD <- read_dta("Data/CONSISTENT_MEDIATION_GROUPC_ODwide.dta")
CONSISTENT_MEDIATION_OD <- read_dta("UKB_paper13_heatmap_updating/CONSISTENT_MEDIATION_GROUPC_ODwide.dta")
dim(CONSISTENT_MEDIATION_OD)
head(CONSISTENT_MEDIATION_OD)
colnames(CONSISTENT_MEDIATION_OD)

CONSISTENT_MEDIATION_OD <- CONSISTENT_MEDIATION_OD %>% 
  dplyr::select(protein, betate, pte, betacde, pcde, betaintref, pintref, betaintmed, pintmed, betapie, ppie, betap_intref, betap_pie) %>% 
  dplyr::rename(Protein = protein, TE = betate, p_TE = pte, CDE = betacde, p_CDE = pcde,
                INTREF = betaintref, p_INT = pintref, INTMED = betaintmed, p_INTMED = pintmed, PIE = betapie, p_PIE = ppie, percent_interaction = betap_intref, percent_mediated = betap_pie) %>% 
  mutate(p_mediated = ifelse(!is.na(p_TE) & p_TE < 0.05, 0.01, 0.1),
         p_interaction = ifelse(!is.na(p_TE) & p_TE < 0.05 & !is.na(p_INT) & p_INT < 0.05, 0.01, 0.1)) %>% 
  mutate_all(~ ifelse(is.na(.) | . == "#DIV/0!" | . == "", NA, .))

CONSISTENT_MEDIATION_OD.long <- to_long(CONSISTENT_MEDIATION_OD, keys = 'term',
                                        values = c('estimate','p'), 
                                        c('TE','CDE','INTREF', 'INTMED', 'PIE', 'percent_interaction', 'percent_mediated'),
                                        c('p_TE','p_CDE','p_INT', 'p_INTMED', 'p_PIE', 'p_interaction', 'p_mediated'))
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ------------------------------------- Part 2 Heat maps --------------------------------------#
# ---------------------------------------------------------------------------------------------#
# ---------------------------------------------------------------------------------------------#
# --------------------------------------
# Heatmap for all outcomes
# --------------------------------------
colnames(NO_MEDIATION_FA)
# combine all outcome

# List of datasets
datasets <- list(NO_MEDIATION_FA.long, INCONSISTENT_MEDIATION_FA.long, CONSISTENT_MEDIATION_FA.long, INCONSISTENT_MEDIATION_OD.long, CONSISTENT_MEDIATION_OD.long)

# List of modality labels
modalities <- c("No mediation (FA)", "Inconsistent mediation (FA)", "Consistent mediation (FA)", "Inconsistent and no mediation combined (OD)", "Consistent mediation (OD)")

# Function to append modality column
append_modality <- function(data, modality) {
  data %>% mutate(p_new = ifelse(p == "<0.001", 0.0009, as.numeric(p)),
                  modality = modality) %>% select(-p)
}

# Apply function to each dataset
datasets_new <- map2(datasets, modalities, append_modality)

# Combine datasets
all_mediation.long <- bind_rows(datasets_new)



heatmap <- function(data, outcome = "No mediation (FA)"){
  # color for mediation
  pal <- scico(7, palette = 'acton')
  
  # assign total number of proteins for different outcomes
  if (outcome %in% c("No mediation (FA)", "Inconsistent mediation (FA)", "Consistent mediation (FA)")){
    number_proteins = 147
  } else if (outcome %in% c("Inconsistent and no mediation combined (OD)", "Consistent mediation (OD)")){
    number_proteins = 147
  }
  
  # set limits based on the maximum absolute value across all exposure-outcome pairs to ensure consistency.
  # and only based on those significant and not mediation
  temp_data <- data %>% filter(term != "percent_mediated" & p_new < 0.05 & !is.na(estimate))
  max_abs_value <- max(abs(c(range(as.numeric(temp_data$estimate)))))
  
  data.long <- data %>% 
    mutate(ap=ifelse(p_new < 0.001, 1,
                     ifelse(p_new >= 0.001 & p_new < 0.01, 2, 
                            ifelse(p_new >= 0.01 & p_new < 0.05, 3, 4))),
           ap=ifelse(is.na(p_new), 4, ap),
           aq=ifelse(p_new < 0.05 , "Pass","insig"),
           aq.mediation=ifelse(p_new < 0.05 , "Pass.mediation", "insig"),
           bg.line=ifelse(term %in% c("percent_mediated", "percent_interaction"), "White", "Dark Grey"),
           bg.color=ifelse(term %in% c("percent_mediated", "percent_interaction"), "Dark Grey", "White"),
           break.mediation = ifelse(!term %in% c("percent_mediated", "percent_interaction"), NA,
                                    ifelse(abs(estimate) <= 0.01, 1,
                                           ifelse(abs(estimate) > 0.01 & abs(estimate) <= 0.05, 2,
                                                  ifelse(abs(estimate) > 0.05 & abs(estimate) <= 0.10, 3,
                                                         ifelse(abs(estimate) > 0.10 & abs(estimate) <= 0.20, 4,
                                                                ifelse(abs(estimate) > 0.20 & abs(estimate) <= 0.50, 5, 6))))))) %>% 
    arrange(factor(Protein, levels = unique(data$Protein)), factor(term, levels = c('TE', 'CDE','INTREF', 'INTMED', 'PIE','percent_interaction','percent_mediated'))) %>% 
    mutate(term = ifelse(term == "percent_mediated", "%mediated", 
                         ifelse(term == "percent_interaction", "%interaction", term)),
           break.mediation = factor(break.mediation, levels = c("1", "2", "3", "4", "5", "6", "NaN"))) %>% 
    filter(modality == outcome) %>% 
    # change the var names to match the footnote
    mutate(term = ifelse(term == "TE", "TERERI", 
                         ifelse(term == "CDE", "ERERI_CDE", 
                                ifelse(term == "INTREF", "ERERI_INTREF", 
                                       ifelse(term == "INTMED", "ERERI_INTMED", 
                                              ifelse(term == "PIE", "ERERI_PIE", term)))))) %>% 
    # remove term == %interaction
    filter(term != "%interaction") 
    
  p.plot <- ggplot(data = data.long, aes(x = forcats::fct_rev(factor(Protein, levels = unique(Protein))), y = factor(term, levels = c('TERERI', 'ERERI_CDE','ERERI_INTREF', 'ERERI_INTMED', 'ERERI_PIE','%mediated'))))+
    geom_tile(color = data.long$bg.line, fill = data.long$bg.color)+
    geom_point(data= subset(data.long, term %in% c('TERERI', 'ERERI_CDE','ERERI_INTREF', 'ERERI_INTMED', 'ERERI_PIE')),
               aes(shape=factor(aq),
                   size=factor(ap), 
                   fill=estimate), na.rm = FALSE)+
    scale_fill_scico(palette = "vik", midpoint = 0, 
                     limits = c(-0.2, max_abs_value+0.05),
                     aesthetics = c("colour","fill")) +
    scale_size_manual(values=c(6, 4, 2.5, 0.5), labels = c("< 0.001", "< 0.01", "< 0.05", "\u2265 .05"))+
    new_scale_color() +
    geom_point(data = subset(data.long, term %in% c("%mediated")), size=6,
               # data = subset(data.long, term %in% c("%interaction", "%mediated")), size=6,
               aes(shape=factor(aq.mediation),
                   color=factor(break.mediation)), na.rm = FALSE)+
    scale_shape_manual(values=c('insig'=1, 'Pass.mediation'=16, 'Pass'=21), drop = FALSE, guide = "none")+
    scale_color_manual(breaks=c(1, 2, 3, 4, 5, 6), drop = FALSE, labels = c("\U2264 1%", "1% - 5%", "5% - 10%", "10% - 20%", "20% - 50%", "> 50%"), values=pal)+
    guides(size = guide_legend(override.aes = list(shape = c(21), fill = c("black")))) +
    labs(title=paste("Heatmap (Four-way decomposition models of selected proteins - ",  outcome, ")", sep = ""),
         # subtitle=paste("", sep = ""),
         x=paste("Proteins (", sum(!is.na(unique(data.long$Protein)))," out of ", number_proteins, " selected proteins)", sep = ""),
         y="",
         size=paste("p-value\nsolid circle: p<.05"), fill=(expression(paste(beta," coefficients"))), color=("% range"),
         # caption="TE: Total effect; CDE: Controlled direct effect; INTREF: Interaction referent;\nINTMED: Mediated interaction; PIE: Pure indirect effect; \n% interaction is the percent of total effect that is pure interaction effect. No p-values were generated; \n% mediated is the percent of total effect that is pure indirect effect. No p-values were generated.") +
         caption="TERERI: Total excess relative risk; ERERI_CDE: Excess relative risk due to neither mediation nor interaction or controlled direct effect;\nERERI_INTREF: Excess relative risk due to interaction only or interaction referent; ERERI_INTMED: Excess relative risk due to mediated interaction or mediated interaction;\nERERI_PIE: Excess relative risk due to mediation only or pure indirect effect; \n% mediated is the percent of total effect that is pure indirect effect. No p-values were generated.") +
    theme(plot.title = element_text(color="Dark blue", size=13, face="bold.italic", hjust = 0.5),
          plot.subtitle=element_text(size=10, hjust=0.5, face="italic", color="Dark blue"),
          plot.caption=element_text(size=9, hjust=0.5, color="Dark grey"),
          axis.title.x = element_text(color="deepskyblue", size=11, face="bold"),
          axis.text.x = element_text(angle = 90, size = 9, vjust = 0.5),
          aspect.ratio=2/7)+
    coord_fixed()
  
  return(p.plot)
}


dir.create(paste(WorkingDirectory,"Output/plot",sep=""), recursive = TRUE)
plot.out.folder <- paste(WorkingDirectory,"Output/plot/",sep="")


outcome.list <- c("No mediation (FA)", "Inconsistent mediation (FA)", "Consistent mediation (FA)", "Inconsistent and no mediation combined (OD)", "Consistent mediation (OD)")
for (outcome in outcome.list){
  if (outcome %in% c("Inconsistent mediation (FA)")){
    width_plot = 11
  } else if (outcome %in% c("Consistent mediation (OD)", "Consistent mediation (FA)")){
    width_plot = 14
  } else if (outcome %in% c("Inconsistent and no mediation combined (OD)")){
    width_plot = 20
  } else if (outcome %in% c("No mediation (FA)")){
    width_plot = 22
  }
  heatmap.hide.unsig <- heatmap(data = all_mediation.long, outcome = outcome)
  ggsave(paste(plot.out.folder, outcome, "_CVH_WMMI_proteins_heatmap_hide_unsig.jpeg",sep=""), heatmap.hide.unsig, width = width_plot, height = 8, units = "in", dpi = 300)
} 