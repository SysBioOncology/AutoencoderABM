path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)

load('data/figure2a.RData')

p1 <- ggplot(figure_tupprol_impkill, aes(x=real, y=estimated)) + 
  geom_abline(slope=1, intercept = 0, linetype = 2, size = 2, color = 'gray') + 
  geom_point(size=3) + 
  xlab('real') + 
  ylab('estimated') + 
  theme_bw() + 
  facet_grid(cols = vars(parameter), scales = 'free') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#a5ba73'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        panel.spacing = unit(1.3, 'lines')) +
  ggpubr::stat_cor(method='pearson', label.y = 1) + 
  ylim(0, 1.1)
p1

p2 <- ggplot(figure_imrwalk, aes(x=real, y=estimated)) + 
  geom_boxplot(fill='#a5ba73') + 
  xlab('real') + 
  ylab('estimated') + 
  theme_bw() + 
  geom_text(data=mean_label, y=0, aes(label=label), size=2.5) + 
  facet_grid(cols = vars(parameter), scales = 'free') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#a5ba73'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold')) + 
  ggpubr:: stat_compare_means(comparisons = my_comparisons,
                              method = 'wilcox',
                              label='p.signif') + 
  ylim(0, 1.4)
p2