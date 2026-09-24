library(ggplot2)
library(data.table)
library(dplyr)
library(ggtext)

#-------------------FST--------------------------

data <- read.table("input.fst", header = FALSE)

#adjust according to the chrmosome names
genome <- data %>% 
  filter(grepl("^SUPER_[1-9]$|^SUPER_1[0-7]$", V1))
genome$V1 <- factor(genome$V1, levels = paste0("SUPER_", 1:17))

col <- as.numeric(sub("SUPER_", "", genome$V1))
super <- genome$V2
max_1 <- c(0)
for (i in 1:17) {
  max_1[i+1] <- max(super[which(col == i)], na.rm = TRUE)
}

len <- numeric(17)
for (i in 1:17) {
  len[i] <- length(super[which(col == i)])
}

sum_1 <- cumsum(max_1[1:17])
sum_super <- rep(sum_1, times = len)
position_total <- super + sum_super


fst <- as.numeric(sub("1:2=", "", genome$V6))

# 5. Criar o Dataframe final limpo (Apenas com os 17 cromossomas)
littorina <- data.frame(
  fst = fst,
  position_total = position_total / 1e6,
  CHROM = genome$V1
)

littorina <- littorina[!is.na(littorina$fst), ]

contig_colors <- rep(c("#3fa15c", "#f74f66"), length.out = 17)

ggplot(littorina, aes(x = position_total, y = fst, color = CHROM)) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_manual(values = contig_colors) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(
    title = "FST values across the genome",
    x = "Genomic Position (Mb)",
    y = "FST") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),    # Remove os traços do eixo X
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_markdown(face = "bold", size = 14, hjust = 0.5), 
    plot.subtitle = element_text(size = 11, hjust = 0.5))

#just for chromosome 11
littorina_chr11 <- data.frame(
  fst = fst[genome$V1 == "SUPER_11"],
  pos_mb = super[genome$V1 == "SUPER_11"] / 1e6
)
littorina_chr11 <- littorina_chr11[!is.na(littorina_chr11$fst), ]

ggplot(littorina_chr11, aes(x = pos_mb, y = fst)) +
  geom_point(color = "#3fa15c", alpha = 0.6, size = 1.5) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(
    title = "FST values - chromosome 11",
    x = "Genomic Position",
    y = "FST") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_markdown(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5),
    axis.text.x = element_text(angle = 0)
  )

#----------------------AF------------------------
data <- read.table("input.af_pwc", header = FALSE)

genome <- data %>% 
  filter(grepl("^SUPER_[1-9]$|^SUPER_1[0-7]$", V1))
genome$V1 <- factor(genome$V1, levels = paste0("SUPER_", 1:17))

col <- as.numeric(sub("SUPER_", "", genome$V1))
super <- genome$V2

max_1 <- c(0)
for (i in 1:17) {
  max_1[i+1] <- max(super[which(col == i)], na.rm = TRUE)
}

len <- numeric(17)
for (i in 1:17) {
  len[i] <- length(super[which(col == i)])
}

sum_1 <- cumsum(max_1[1:17])
sum_super <- rep(sum_1, times = len)
position_total <- super + sum_super

freq <- as.numeric(genome$V9)

littorina <- data.frame(
  freq = freq,
  position_total = position_total / 1e6,
  CHROM = genome$V1
)

littorina <- littorina[!is.na(littorina$freq), ]

contig_colors <- rep(c("#3fa15c", "#f74f66"), length.out = 17)

ggplot(littorina, aes(x = position_total, y = freq, color = CHROM)) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_manual(values = contig_colors) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Allele frequency differencies across the genome",
    x = "Genomic Position (Mb)",
    y = "AF") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_markdown(face = "bold", size = 22, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

littorina_chr11 <- data.frame(
  freq = freq[genome$V1 == "SUPER_11"],
  pos_mb = super[genome$V1 == "SUPER_11"] / 1e6)

littorina_chr11 <- littorina_chr11[!is.na(littorina_chr11$freq), ]

ggplot(littorina_chr11, aes(x = pos_mb, y = freq)) +
  geom_point(color = "#3fa15c", alpha = 0.6, size = 1.5) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = "Allele frequency differencies - chromosome 11",
    x = "Genomic Position",
    y = "AF") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_markdown(face = "bold", size = 22, hjust = 0.5),
    axis.text.x = element_text(angle = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

#------------------pi and Tajima'D-----------------------
data <- read.table("input.D", header = FALSE)
#REPEAT FOR PI

genome <- data %>% 
  filter(grepl("^SUPER_[1-9]$|^SUPER_1[0-7]$", V1))
genome$V1 <- factor(genome$V1, levels = paste0("SUPER_", 1:17))

col <- as.numeric(sub("SUPER_", "", genome$V1))

super <- genome$V2

max_1 <- c(0)
for (i in 1:17) {
  max_1[i+1] <- max(super[which(col == i)], na.rm = TRUE)
}

len <- numeric(17)
for (i in 1:17) {
  len[i] <- length(super[which(col == i)])
}

sum_1 <- cumsum(max_1[1:17])
sum_super <- rep(sum_1, times = len)
position_total <- super + sum_super
D <- as.numeric(genome$V5)

littorina <- data.frame(
  D = D,
  position_total = position_total / 1e6,
  CHROM = genome$V1
)
littorina <- littorina[!is.na(littorina$D), ]

contig_colors <- rep(c("#3fa15c", "#f74f66"), length.out = 17)

ggplot(littorina, aes(x = position_total, y = D, color = CHROM)) +
  geom_point(alpha = 0.5, size = 1) + 
  scale_color_manual(values = contig_colors) +
  coord_cartesian(ylim = c(-2.5, 2)) +
  labs(
    title = "Tajima's D values across the genome"
    x = "Genomic Position (Mb)",
    y = "Tajima's D") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 0),
    plot.title = element_markdown(face = "bold", size = 22, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

littorina_chr11 <- data.frame(
  D = D[genome$V1 == "SUPER_11"],
  pos_mb = super[genome$V1 == "SUPER_11"] / 1e6
)

littorina_chr11 <- littorina_chr11[!is.na(littorina_chr11$D), ]

ggplot(littorina_chr11, aes(x = pos_mb, y = D)) +
  geom_point(color = "#3fa15c", alpha = 0.6, size = 1.5) +
  coord_cartesian(ylim = c(-2.5, 2)) +
  labs(
    title = "Tajima's D values - chromosome 11",
    x = "Genomic Position",
    y = "Tajima's D") +
  theme_minimal() +
  theme(
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_markdown(face = "bold", size = 22, hjust = 0.5),
    axis.text.x = element_text(angle = 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12))

#-----------------Wilcoxon--------------------

#REPEAT THIS FOR ALL STATITICS, FOR ALL SPECIES, with different recognizable names
#ADJUST THE INITIAL AND END POSITIONS WITH THE POSITIONS WITH HIGH FST
data_11 <- data[which(data$V1 == "SUPER_11"), ]
data_wilc <- data_11 %>% mutate(region = if_else((V2 >= 33825000 & V2 <= 35750000) | (V2 >= 41887500 & V2 <= 42400000),"inside", "outside"))

#add here all matrix
all_dataframes <- list(
  data_wilc = data_wilc
  ...
  ...
)

wilcoxon <- function(df_name, df) {
  #here, fst values were on the 9th column, af values were on the 6th and pi and tajima's D were on the 5th column.
  column <- if (grepl("af", df_name)) {
    "V9"
  } else if (grepl("fst", df_name)) {
    "V6"
  } else {
    "V5"
  }

  df_clear <- df %>% 
    mutate(across(all_of(column), as.numeric)) %>% 
    filter(!is.na(.[[column]]), !is.na(region))
  
  outside <- df_clear %>% filter(region == "outside") %>% pull(!!sym(column))
  inside   <- df_clear %>% filter(region == "inside") %>% pull(!!sym(column))
  
  test <- wilcox.test(, fora)
  
  return(data.frame(
    Dataframe = df_name,
    Col = column,
    N_in = length(inside),
    N_out = length(outside),
    W = test$statistic,
    p_value = test$p.value,
    Message = "OK"
  ))
}

results <- do.call(rbind, lapply(names(all_dataframes), function(name) {
  rodar_wilcoxon(name, all_dataframes[[name]])
}))

print(results, row.names = FALSE)

