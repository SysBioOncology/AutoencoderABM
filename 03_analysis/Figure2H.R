path <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(path)

library(ggplot2)
library(ggthemes)
library(patchwork)
library(dplyr)

load('data/figure2h.RData')
load('data/rnaseq_data.RData') # can be used to further investigate gene expression

# trained on tcga autoencoder 
stat_comparison <- ggplot(plot_means, aes(x=gene, y=mean, color=category)) + 
  facet_wrap(~parameter, scales='free') + 
  geom_point(size=3, position=position_dodge(width=0.7)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.7),
                alpha=0.5) + 
  scale_color_manual(values=c('#eab139', '#d36f4b', '#843a20')) +
  theme_bw() + 
  theme(strip.background = element_rect(fill='#f3d291'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', hjust=0.5),
        legend.position = 'bottom',
        axis.title.x = element_text(margin=margin(10,0,0,0))) + 
  ylab('normalized\ngene expression') + 
  ggpubr::geom_pwc(data = sign_test_tcga,
                   aes(group = category, y=expression, x = gene), tip.length = 0,
                   method = "wilcox_test", label = "p.signif", vjust=0.1,
                   p.adjust.method='BH',inherit.aes = F, y.position = 1.1,
                   label.size=3 
  ) +
  ylim(-1.75, 1.8) +
  ggtitle('Autoencoder trained on TCGA')
stat_comparison

# trained on synthetic data 
stat_comparison_synthetic <- ggplot(plot_means, aes(x=gene, y=mean, color=category)) + 
  facet_wrap(~parameter, scales='free') + 
  geom_point(size=3, position=position_dodge(width=0.7)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), position=position_dodge(width=0.7),
                alpha=0.5) + 
  scale_color_manual(values=c('#eab139', '#d36f4b', '#843a20')) +
  theme_bw() + 
  theme(strip.background = element_rect(fill='#f3d291'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold', hjust=0.5),
        legend.position = 'bottom',
        axis.title.x = element_text(margin=margin(10,0,0,0))) + 
  ylab('normalized\ngene expression') + 
  ggpubr::geom_pwc(data = sign_test_synthetic,
                   aes(group = category, y=expression, x = gene), tip.length = 0,
                   method = "wilcox_test", label = "p.signif", vjust=0.1,
                   p.adjust.method='BH',inherit.aes = F, y.position = 1.1,
                   label.size=3
  ) +
  ylim(-1.75, 1.8) + 
  ggtitle('Autoencoder trained on synthetic data')
stat_comparison_synthetic

# histogram of parameter values 
hist_plot <- ggplot(all_hist, aes(x=estimated, fill=category)) + 
  geom_histogram(alpha=0.9) + 
  facet_wrap(~parameter) + 
  theme_bw() + 
  theme(strip.background = element_rect(fill='#f3d291'),
        strip.text = element_text(face='bold', size=10),
        axis.title = element_text(face='bold', size=10),
        plot.title = element_text(face='bold'),
        legend.position = 'None') + 
  scale_fill_manual(values=c('#eab139', '#d36f4b', '#843a20')) +
  scale_y_continuous(breaks = c(0, 250, 500))
hist_plot

hist_plot / (stat_comparison + theme(legend.position = 'None')) / stat_comparison_synthetic + 
  plot_layout(heights = c(1, 4, 4))
ggsave(file=paste('output/compilation_gene_expression_comparisons.png', sep=''), 
       dpi=600, width=7, height=9)