path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)
library(readxl)
library(ggh4x)

load('data/figure2beg.RData')


# COMPLEXITY --------------------------------------------------------
complexity_plot <- ggplot(subset(complexity_all, category == 'synthetic'), aes(x=`complexity (training)`, y=`complexity (simulation)`,
                                                                               color=category, fill = category)) + 
  geom_point(size=4) + 
  geom_smooth(method = 'lm', alpha=0.2, linetype = 2) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        legend.position = 'None',
        panel.spacing = unit(1, 'lines')) +
  ggpubr::stat_cor(method='pearson', color = 'black')  +
  scale_color_manual(values=c('#a5ba73', '#046c9a', '#eab139')) + 
  scale_fill_manual(values=c('#a5ba73', '#046c9a', '#eab139')) + 
  xlab('complexity (original)')
complexity_plot 
ggsave(complexity_plot, file='output/optimized_complexity_synthetic.png', 
       dpi=600, width=2.5, height=2.5)

complexity_plot <- ggplot(subset(complexity_all, category == 'tumoroid'), aes(x=`complexity (training)`, y=`complexity (simulation)`,
                                                                              color=category, fill = category)) + 
  geom_point(size=4) + 
  geom_smooth(method = 'lm', alpha=0.2, linetype = 2) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        legend.position = 'None',
        panel.spacing = unit(1, 'lines')) +
  ggpubr::stat_cor(method='pearson', color = 'black')  +
  scale_color_manual(values=c('#046c9a', '#eab139')) + 
  scale_fill_manual(values=c('#046c9a', '#eab139')) + 
  xlab('complexity (original)')
complexity_plot 
ggsave(complexity_plot, file='output/optimized_complexity_tumoroid.png', 
       dpi=600, width=2.5, height=2.5)


# NEIGHBORS ------------------------------------------------------------------
neighbors_plot <- ggplot(subset(neighbors_all, category == 'synthetic'), 
                         aes(x=`mean (training)`, y=`mean (simulation)`, color=category,
                                                                             fill = category)) + 
  geom_point(size=3, alpha=0.6) + 
  geom_smooth(method = 'lm', alpha=0.2, linetype = 2) + 
  theme_bw() + 
  facet_grid2(rows = vars(`neighbor_cell_type`),
              scales = 'free', independent = 'all') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', hjust=0.5),
        panel.spacing = unit(1, 'lines')) +
  ggpubr::stat_cor(method='pearson', color='black') + 
  scale_color_manual(values=c('#a5ba73')) + 
  scale_fill_manual(values=c('#a5ba73')) + 
  xlab('mean count (original)') + 
  ylab('mean count (simulation)')
neighbors_plot
ggsave(neighbors_plot, file='output/optimized_neighbors_synthetic.png', 
       dpi=600, width=3, height=4)


neighbors_plot <- ggplot(subset(neighbors_all, category == 'tumoroid'), aes(x=`mean (training)`, y=`mean (simulation)`, color=category,
                                                                            fill = category)) + 
  geom_point(size=3, alpha=0.6) + 
  geom_smooth(method = 'lm', alpha=0.2, linetype = 2) + 
  theme_bw() + 
  facet_grid2(rows = vars(`neighbor_cell_type`),
              scales = 'free', independent = 'all') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', hjust=0.5),
        panel.spacing = unit(1, 'lines')) +
  ggpubr::stat_cor(method='pearson', color='black') + 
  scale_color_manual(values=c('#046c9a')) + 
  scale_fill_manual(values=c('#046c9a')) + 
  xlab('mean count (original)') + 
  ylab('mean count (simulation)')
neighbors_plot
ggsave(neighbors_plot, file='output/optimized_neighbors_tumoroid.png', 
       dpi=600, width=3, height=4)




neighbors_plot <- ggplot(subset(neighbors_all, category == 'TCGA'), aes(x=`mean (training)`, y=`mean (simulation)`, color=category,
                                                                        fill = category)) + 
  geom_point(size=3, alpha=0.6) + 
  geom_smooth(method = 'lm', alpha=0.2, linetype = 2) + 
  theme_bw() + 
  facet_grid2(rows = vars(`neighbor_cell_type`),
              scales = 'free', independent = 'all') + 
  theme(legend.position = 'None',
        strip.background = element_rect(fill='white'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', hjust=0.5),
        panel.spacing = unit(1, 'lines')) +
  ggpubr::stat_cor(method='pearson', color='black') + 
  scale_color_manual(values=c('#eab139')) + 
  scale_fill_manual(values=c('#eab139')) + 
  xlab('mean count (original)') + 
  ylab('mean count (simulation)')
neighbors_plot
ggsave(neighbors_plot, file='output/optimized_neighbors_tcga.png', 
       dpi=600, width=3, height=4)