library("tidyverse")
library("htmlwidgets")
library("manhattanly")
GWAS_result1 <- read.csv("C://Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/impute_out/Mg/output.csv")
GWAS_result1[GWAS_result1==Inf] <- NA
GWAS_result1[GWAS_result1==-Inf] <- NA
GWAS_result<-GWAS_result1[complete.cases(GWAS_result1),]
# Calculate the BHY threshold
m <- nrow(GWAS_result)
GWAS_result <- GWAS_result[order(GWAS_result$absolute_theta),]
s <- 1.0
i <- 0
for (p in GWAS_result$absolute_theta) {
  p
  i <- i+1
  if (i > 1) {
    s <- s + 1.0/(i-1)
  }
  thes_pval <- ((i + 1.0) / m) * 0.05 / s
  if (p > thes_pval) {break
  }
}
thes_pval_original <- thes_pval
bhy_thres <- -log10(thes_pval)
# calculate bonferroni_threshold
bt <- 0.05 / (nrow(GWAS_result)*1135) # times max number of tests per p-value
bf_thres <- -log10(bt)
data_cum <- GWAS_result %>% 
  group_by(CHROM) %>% 
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHROM, bp_add)
GWAS_result <- GWAS_result %>% 
  inner_join(data_cum, by = "CHROM") %>% 
  mutate(bp_cum = POS + bp_add)
axis_set <- GWAS_result %>% 
  group_by(CHROM) %>% 
  summarize(center = mean(bp_cum))
ylim <- abs(floor(log10(min(GWAS_result$absolute_theta)))) +1
png("C://Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/impute_out/Mg/output_absolute_theta.png", bg = "white", width = 9.75, height = 3.25, units = "in", res = 1200, pointsize = 4)
manhplot <- ggplot(GWAS_result, aes(x = bp_cum, y = (absolute_theta), # size = 1, 
                                  color = as_factor(CHROM))) +
  geom_point(alpha = 0.5) +
geom_hline(yintercept = bf_thres, color = "red", linetype = "dashed") +
geom_hline(yintercept = bhy_thres, color = "blue", linetype = "dashed") +
  scale_x_continuous(label = axis_set$CHROM, breaks = axis_set$center) +
  labs(x = expression(paste("SNP Position on ", italic("C. excelsa"), " Reference Genome")), 
       y = "Absolute Theta") + 
  theme_minimal() +
  guides(colour="none")
  theme(
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
  )
print(manhplot)
dev.off()
attach(GWAS_result)
GWAS_result <- GWAS_result[order(-absolute_theta),]
detach(GWAS_result)
top_perc <- (nrow(GWAS_result)/100) * 0.01
top_perc <-round(top_perc, digits = 0)
GWAS_result <- head(GWAS_result,top_perc)
names(GWAS_result)[names(GWAS_result) == "CHROM"] <- "CHR"
GWAS_result$P <- (0.0000001/GWAS_result$absolute_theta)
names(GWAS_result)[names(GWAS_result) == "POS"] <- "BP"
GWAS_result$CHR <- gsub("[^0-9]", "", GWAS_result$CHR)
GWAS_result$CHR <- as.numeric(GWAS_result$CHR)


html_file <- manhattanly(GWAS_result, annotation1 = "CHR", annotation2 = "BP")
saveWidget(html_file, "C://Users/Oliver/Desktop/UoN/Individual_project/GWAS_out/impute_out/Mg/output_absolute_theta.html", selfcontained = T)
