path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)

load('data/suppfig3c.RData')
load('data/suppfig3-2.RData')

ggplot(plot_comparison, aes(x=`cell type`, y=count, fill=`immune subtype`)) + 
  geom_boxplot(width=0.5) + 
  ggpubr::stat_compare_means(aes(group = `immune subtype`), label = "p.signif",
                             label.y.npc = 0.97) + 
  theme_bw() + 
  theme(legend.position = 'bottom', 
        axis.title = element_text(face='bold', size=12)) + 
  scale_fill_manual(values = c('#eab139', '#d36f4b'))
ggsave('output/tcga_cell_type_per_immune_subtype.png', dpi=600, 
       height=3.5, width=3.5)

ggplot(comp_par, aes(x=TCGA, y=synthetic)) + 
  geom_smooth(method='lm', color = 'gray', fill = 'gray') + 
  geom_point() + 
  facet_wrap(~parameter) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill='#f3d291'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        legend.position = 'None') 
ggsave('output/tcga_comparison_parameters.png', dpi=600, height=2.5, width=8)