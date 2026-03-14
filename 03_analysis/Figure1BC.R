path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(readxl)
library(ggh4x)

# COMPLEXITY -------------------------------
# Get data ---------------------------------------------------------------------
folder = 'data/tumoroid_'
complexity <- read.csv(paste(folder, 'complexity_scores_validation.csv', sep=''), check.names=F)
complexity$data <- 'tumoroid'

folder = 'data/synthetic_'
temp <- read.csv(paste(folder, 'complexity_scores_validation.csv', sep=''), check.names=F)
temp$data <- 'synthetic'
temp <- temp[!round(temp$`complexity (training)`, 2) %in% c(1.27),]

complexity <- rbind(complexity, temp)
complexity <- complexity[!is.na(complexity$`complexity (reconstruction)`),]


# Create some plots ------------------------------------------------------------
complexity$data <- factor(complexity$data, levels=c('synthetic', 'tumoroid', 'TCGA'))
xy_line <- data.frame('x'=c(1, 2), 'y'=c(1, 2))

complexity <- ggplot(subset(complexity, `complexity (training)` < 2 & `complexity (reconstruction)` < 2), 
                     aes(x=`complexity (training)`, y=`complexity (reconstruction)`, color=data)) + 
  geom_line(data=xy_line, aes(x=x, y=y), inherit.aes=F,
            linetype=2, color='darkgray', size=1) + 
  geom_point(size=3, alpha=0.4) + 
  facet_grid2(rows = vars(data), scales = 'free',
              independent = 'all') + 
  theme_bw() + 
  theme(strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=11),
        axis.title = element_text(face='bold', size=11),
        plot.title = element_text(face='bold'),
        legend.position = 'None',
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.5, "lines")) +
  ggpubr::stat_cor(method='pearson', color='black', size=3) + 
  xlab('complexity (original)') + 
  ylab('complexity (reconstruction)') + 
  xlim(1, 2) + 
  scale_color_manual(values=c('#a5ba73', '#046c9a', '#eab139'))
complexity

ggsave(file=paste('output/complexity_correlation.png', sep=''), dpi=600, height=5, width=3)


# NEIGHBORS -------------------------------
# Get data ---------------------------------------------------------------------
folder = 'data/tcga_'
neighbors_mean <- read.csv(paste(folder, 'neighbors_mean_validation.csv', sep=''),
                           check.names=F)
neighbors_mean$data <- 'TCGA'

folder = 'data/tumoroid_'
temp <- read.csv(paste(folder, 'neighbors_mean_validation.csv', sep=''),
                 check.names=F)
temp$data <- 'tumoroid'
neighbors_mean <- rbind(neighbors_mean, temp)


folder = 'data/synthetic_'
temp <- read.csv(paste(folder, 'neighbors_mean_validation.csv', sep=''),
                 check.names=F)
temp$data <- 'synthetic'
neighbors_mean <- rbind(neighbors_mean, temp)


# Create some plots ------------------------------------------------------------
neighbors_mean <- neighbors_mean[!is.na(neighbors_mean$`mean (reconstruction)`),]
neighbors_mean$data <- factor(neighbors_mean$data, 
                              levels=c('synthetic', 'tumoroid', 'TCGA'))
neighbors <- ggplot(neighbors_mean, aes(x=`mean (training)`, y=`mean (reconstruction)`,
                                        color = data)) + 
  geom_point(size=3, alpha=0.4) + 
  facet_grid2(cols = vars(neighbor_cell_type), 
              rows = vars(data), scales = 'free',
              independent = 'all') + 
  theme_bw() + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=11),
        axis.title = element_text(face='bold', size=11),
        plot.title = element_text(face='bold', hjust=0.5),
        panel.spacing.x = unit(0.5, "lines"),
        panel.spacing.y = unit(0.5, "lines")) +
  #ggpubr::stat_cor(method='spearman', label.y.npc = 0.5, label.x.npc = 0.5) + 
  ggpubr::stat_cor(method='pearson', color='black', size=3) + 
  ggtitle('Neighbors of tumor cells') + 
  scale_color_manual(values=c('#a5ba73', '#046c9a', '#eab139')) + 
  xlab('mean count (original)') + 
  ylab('mean count (simulation)')
neighbors
ggsave(file=paste('output/neighbors_correlation.png', sep=''), dpi=600, height=6, width=5)

