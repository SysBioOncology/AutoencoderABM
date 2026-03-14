path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(readxl)

load('data/figure2d.RData')

estimated_plot <- ggplot(subset(all_results_both, `cell type` == 'tumor cells'), 
                         aes(y=IMpkill, x=experiment, fill = experiment)) + 
  geom_boxplot(width=0.5, outlier.shape = NA) + 
  geom_jitter(size=2, width=0.2) + 
  scale_color_manual(values = c('#eccbae', '#ABDDDE')) + 
  scale_fill_manual(values = c('#eccbae', '#ABDDDE')) + 
  theme_bw() + 
  ylab('estimated IMpkill') + 
  ggpubr::geom_pwc(
    aes(group = experiment), tip.length = 0,
    method = "wilcox_test", label = "p.signif",
    p.adjust.method='BH',
    label.size=4
  ) + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#046C9A'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        axis.text.x = element_text(angle=15, vjust=1, hjust=1)) +
  ylim(0, 0.8)
estimated_plot

ggsave(file='output/tumoroid_comparison_optimized_parameters.png', dpi=600, height=3, width=3) 

