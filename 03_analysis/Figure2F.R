path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(readxl)

load('data/figure2f.RData')

ggplot(plot_results_long, aes(x=autoencoder, y=estimated, fill=`immune subtype`)) + 
  geom_boxplot(width=0.7) + 
  facet_wrap(~parameter) + 
  scale_fill_manual(values = c('#eab139', '#d36f4b')) + 
  theme_bw() + 
  ylab('estimated') + 
  theme(legend.position = 'bottom',
        strip.background = element_rect(fill='#f3d291'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(hjust=0.5, size=12, face='bold'),
        panel.spacing = unit(1, 'lines')) + 
  xlab('data used for autoencoder optimization') +
  ylab('estimated') + 
  ggpubr::geom_pwc(
    aes(group = `immune subtype`), tip.length = 0,
    method = "wilcox_test", label = "p.signif",
    p.adjust.method='BH',
    label.size=4
  ) +
  ylim(0, 1.15)
ggsave('output/tcga_comparison_autoencoders.png', dpi=600, width=5, height=3)