
# 4.Stats #

#### Package setup ####

library(openxlsx)
library(dplyr)
library(corrplot)
library(Hmisc)
library(ggplot2)
library(ggstatsplot)
library(ggside)

theme_set(theme_bw())

#### End ####

#### Input and reformat data ####
corr_data = read.xlsx("Results/4.Stats/correlations_tab.xlsx", rowNames = T, sheet = "Final_corr_matrix")

# Select only numeric columns 
my_data =  corr_data %>% select_if(is.numeric)
my_data = my_data[,colnames(my_data) != "BDG_date"]

colnames(corr_data)[!colnames(corr_data) %in% colnames(my_data)] #Deleted columns

# Reformat names for plot
colnames(my_data) <- gsub("\\.|_", " ", colnames(my_data))

#### End ####

#### Check normality with Shapiro test ####

shapiro <- list()
for(i in colnames(my_data)){
  shapiro[[i]] <- cbind("variable" = i,
                      
                      "shapiro.W.Stat" = shapiro.test(my_data[,i])$statistic,
                      "shapiro.p.value" = shapiro.test(my_data[,i])$p.value)
  
}
shapiro <- as.data.frame(do.call(rbind, shapiro))

# All variables are nonparametric except for age
#### End ####

#### Spearman correlation ####

res <- rcorr(as.matrix(my_data), type = "spearman")

# Extract the correlation coefficients
#res$r
# Extract p-values
#res$P

# Plot correlations
col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))

corrplot(res$r, method = "circle", col = rev(col(200)),  
         type = "lower", #order = "hclust", 
         # addCoef.col = "black", number.cex = 0.7, number.font = 1,# Add coefficient of correlation
         tl.col = "black", tl.srt = 90, tl.cex = 0.8,#Text label color and rotation
         # Combine with significance level
         p.mat = res$P, sig.level = 0.05, insig = "blank", na.label.col = "white",
         # hide correlation coefficient on the principal diagonal
         diag = FALSE 
)

#### End ####

#### Single correlations ####

## Fungal ITS2 reads - 6th day
ggscatterstats(
  data = my_data,
  x = "Fungal ASV counts", 
  y = '6d evaluation', 
  xlab = "Fungal ITS2 read count",
  ylab = "6th day evaluation",
  point.label.args = list(alpha = 0.7, size = 4, color = "grey50"),
  xfill = "#CC79A7", ## fill for marginals on the x-axis
  yfill = "#009E73", ## fill for marginals on the y-axis
  type ="n",
#  results.subtitle = F
) + theme(text = element_text(size = 15),
          plot.title = element_text(size = 18))

#### End ####

#### Mann-Whitney ####

mwut <- list()

for (n in c("Fungal ASV counts")){
  for (i in c("ICU discharge (alive/deceased 1/0)", "Intestinal ischemia", "Tumor", "Peritonitis", "Cirrhosis", "Antibiotic therapy (+5d)", "Antifungal (+5d)", "Antibiotic therapy (-14d)", "Bacterial growth", "Fungal growth")){
    U_test <- wilcox.test(my_data[,n] ~ my_data[,i], data = my_data, paired = F)
    z <- abs(qnorm(U_test$p.value/2))
    tab <- c(U_test$method, paste("data:", n, "and", i), U_test$statistic, U_test$p.value, z)
    mwut[[paste0(i, "_", n)]] <- tab
  }
}

mwut <- as.data.frame(do.call(rbind, mwut))
colnames(mwut) <- c("Test", "Variables", "Statistic", "p-value", "z")

# Plot significant ones

# Peritonitis
my_data <- my_data %>%
  mutate(Peritonitis_character = ifelse(Peritonitis == 0, "No", "Yes"))

ggbetweenstats(
  data  = my_data,
  x     = "Peritonitis_character",
  y     = 'Fungal ASV counts',
  type = "n",
  p.adjust.method = "BH",
  results.subtitle = F,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.4, size = 6, stroke = 0, na.rm = TRUE),
  centrality.plotting = F
) + xlab("Peritonitis") + ylab("Fungal ITS2 read count") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 18))

# ICU discharge
my_data <- my_data %>%
  mutate(`ICU_discharge` = ifelse(`ICU discharge (alive/deceased 1/0)` == 0, "Deceased", "Alive"))

ggbetweenstats(
  data  = my_data,
  x     = 'ICU_discharge',
  y     = 'Fungal ASV counts',
  type = "n",
  p.adjust.method = "BH",
  results.subtitle = F,
  point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.6), alpha =
                      0.4, size = 6, stroke = 0, na.rm = TRUE),
  centrality.plotting = F
) + xlab("Discharge from ICU") + ylab("Fungal ITS2 read count") +
  theme(text = element_text(size = 15),
        plot.title = element_text(size = 18))

#### End ####

#### Metadata metrics ####

# Calculate mean, median, maximum and minimum
metrics = list()
for(v in colnames(my_data)){
  if(is.numeric(my_data[,v])){
    mean = mean(my_data[,v], na.rm = T)
    median = median(my_data[,v], na.rm = T)
    max = max(my_data[,v], na.rm = T)
    min = min(my_data[,v], na.rm = T)
    metrics[[v]] = c(mean, median, max, min)
    names(metrics[[v]]) = c("mean", "median", "max", "min")
  }}
metrics <- as.data.frame(do.call(rbind, metrics)) 

# Calculate metrics per group of culture and ITS2 positivity

# Recover character variables
my_data <- corr_data

my_data$category_fungi = as.factor(my_data$category_fungi)

metrics_grouped = list()
for(v in colnames(my_data)){
  if(is.numeric(my_data[,v])){
    Q = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = quantile, na.rm = T)
    IQR = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = IQR, na.rm = T)[,"x"]
    mean = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = mean, na.rm = T)[,"x"]
    median = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = median, na.rm = T)[,"x"]
    min = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = min, na.rm = T)[,"x"]
    max = aggregate(x= my_data[,v], by = list(my_data$category_fungi), FUN = max, na.rm = T)[,"x"]
    metrics_grouped[[v]] <- cbind(v, Q, IQR, mean, median, min, max)
    names(metrics_grouped[[v]]) <- c("v", "category", "Percentile", "IQR", "mean", "median", "min", "max")
  }}
metrics_grouped <- as.data.frame(do.call(rbind, metrics_grouped))

# Calculate number and percentage for character variables per group of culture and ITS2 positivity

sex <- my_data %>% group_by(category_fungi, sex) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(sex, count, fill = 0) 
sex$percentage <- sex[,"m"] / (sex[,"m"] + sex[,"w"]) * 100
Alcoholism <- my_data %>% group_by(category_fungi, `Alcoholism`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`Alcoholism`, count, fill = 0) 
Alcoholism$percentage <- Alcoholism[,"1"] / (Alcoholism[,"1"] + Alcoholism[,"0"]) * 100
Smoking <- my_data %>% group_by(category_fungi, Smoking) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Smoking, count, fill = 0) 
Smoking$percentage <- Smoking[,"1"] / (Smoking[,"1"] + Smoking[,"0"]) * 100
`ICU.discharge.(alive/dead.1/0)` <- my_data %>% group_by(category_fungi, `ICU.discharge.(alive/dead.1/0)`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`ICU.discharge.(alive/dead.1/0)`, count, fill = 0) 
`ICU.discharge.(alive/dead.1/0)`$percentage <- `ICU.discharge.(alive/dead.1/0)`[,"1"] / (`ICU.discharge.(alive/dead.1/0)`[,"1"] + `ICU.discharge.(alive/dead.1/0)`[,"0"]) * 100
Intestinal_ischemia <- my_data %>% group_by(category_fungi, Intestinal_ischemia) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Intestinal_ischemia, count, fill = 0) 
Intestinal_ischemia$percentage <- Intestinal_ischemia[,"1"] / (Intestinal_ischemia[,"1"] + Intestinal_ischemia[,"0"]) * 100
Tumor <- my_data %>% group_by(category_fungi, Tumor) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Tumor, count, fill = 0) 
Tumor$percentage <- Tumor[,"1"] / (Tumor[,"1"] + Tumor[,"0"]) * 100
Peritonitis <- my_data %>% group_by(category_fungi, Peritonitis) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Peritonitis, count, fill = 0) 
Peritonitis$percentage <- Peritonitis[,"1"] / (Peritonitis[,"1"] + Peritonitis[,"0"]) * 100
Cirrhosis <- my_data %>% group_by(category_fungi, Cirrhosis) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Cirrhosis, count, fill = 0) 
Cirrhosis$percentage <- Cirrhosis[,"1"] / (Cirrhosis[,"1"] + Cirrhosis[,"0"]) * 100
`Antibiotic.therapy.(+5d)` <- my_data %>% group_by(category_fungi, `Antibiotic.therapy.(+5d)`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`Antibiotic.therapy.(+5d)`, count, fill = 0) 
`Antibiotic.therapy.(+5d)`$percentage <- `Antibiotic.therapy.(+5d)`[,"1"] / (`Antibiotic.therapy.(+5d)`[,"1"] + `Antibiotic.therapy.(+5d)`[,"0"]) * 100
`Antibiotic.therapy.(-14d)` <- my_data %>% group_by(category_fungi, `Antibiotic.therapy.(-14d)`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`Antibiotic.therapy.(-14d)`, count, fill = 0) 
`Antibiotic.therapy.(-14d)`$percentage <- `Antibiotic.therapy.(-14d)`[,"1"] / (`Antibiotic.therapy.(-14d)`[,"1"] + `Antibiotic.therapy.(-14d)`[,"0"]) * 100
`Antifungal.(+5d)` <- my_data %>% group_by(category_fungi, `Antifungal.(+5d)`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`Antifungal.(+5d)`, count, fill = 0) 
`Antifungal.(+5d)`$percentage <- `Antifungal.(+5d)`[,"1"] / (`Antifungal.(+5d)`[,"1"] + `Antifungal.(+5d)`[,"0"]) * 100
`Bloodculture.(+-5d)` <- my_data %>% group_by(category_fungi, `Bloodculture.(+-5d)`) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(`Bloodculture.(+-5d)`, count, fill = 0) 
`Bloodculture.(+-5d)`$percentage <- `Bloodculture.(+-5d)`[,"1"] / (`Bloodculture.(+-5d)`[,"1"] + `Bloodculture.(+-5d)`[,"0"]) * 100
Bacterial_growth <- my_data %>% group_by(category_fungi, Bacterial_growth) %>% summarise(count = n(), .groups = 'drop') %>%
  tidyr::spread(Bacterial_growth, count, fill = 0) 
Bacterial_growth$percentage <- Bacterial_growth[,"1"] / (Bacterial_growth[,"1"] + Bacterial_growth[,"0"]) * 100

#### End ####
