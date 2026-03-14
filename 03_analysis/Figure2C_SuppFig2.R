path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(readxl)

# get data for figure 2C and supplementary figure 2
load('data/figure2c_supp.RData')

# Figure 2C -------------------------------------------------------------------
p1 <- ggplot(best_results_long, aes(x=real, y=estimated)) + 
  geom_abline(slope=1, intercept = 0, linetype = 2, size = 1, color = 'gray') + 
  geom_point(size=2, show.legend=F, alpha=0.6) + 
  xlab('real') + 
  ylab('estimated') + 
  theme_bw() + 
  facet_grid(cols = vars(parameter), scales = 'free') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#a5ba73'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        panel.spacing = unit(1.3, 'lines'),
        axis.title.x = element_blank()) +
  ggpubr::stat_cor(method='pearson', label.y = 1.1, vjust=2)+ 
  ylim(0, 1.1) 
p1

p3 <- ggplot(best_results_long, aes(x=`category (real)`, y=estimated)) + 
  geom_boxplot(width=0.2, fill = '#7f9549') +
  theme_bw() + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#a5ba73'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        axis.text.x = element_text(angle=15, vjust=1, hjust=1),
        panel.spacing = unit(1.3, 'lines')) + 
  xlab('real')
p3

p1 / p3 + plot_layout(heights = c(2, 1))
ggsave(file=paste('output/compilation_fixed_tupprol.png', sep=''), 
       dpi=600, height=3, width=3)


# Supplementary figure 2 ------------------------------------------------------
inv <- ggplot(investigate, aes(x=`TUpprol (real)`, y=diff)) + 
  geom_hline(yintercept=0, color='black', linetype=2, size=1) + 
  geom_violin(fill = '#e0e7cf') + 
  geom_boxplot(width=0.1, fill = '#7f9549') + 
  theme_bw() + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='#a5ba73'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10)) + 
  ylab('real-estimated\nIMpkill') + 
  coord_flip() 
inv
ggsave(file=paste('output/supplementary_figure2.png', sep=''), 
       dpi=600, height=8, width=2)
