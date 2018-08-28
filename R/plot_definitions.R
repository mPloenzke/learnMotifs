#' Plot effect size versus activation difference.
#'
#' Plots the estimated motif effect size (beta) against the mean between-class activation difference.
#'
#' @param betas Tibble containing effect size information, filter name, mean activation difference, and information.
#' @param combo Whether both de novo and annotated motifs are included.
#' @param fl Filename to save plots to. If not provided, a plot object will be returned.
#'
#' @return A ggplot object.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords plot beta betas
#'
#' @examples
#' a <- runif(10,0,1)
#' b <- rnorm(10)
#' c <- paste(1:10)
#' d <- rnorm(10,1)
#' mat <- dplyr::tibble(Difference=a,Effect=b,Filter=c,Information=d)
#' plot_activation_difference(mat)
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_hline geom_vline scale_size theme ylab xlab ggsave element_rect element_text
#' @importFrom ggrepel geom_text_repel
#' @export
plot_activation_difference <- function(betas,combo=F,fl=NULL) {
  if (!combo) {
    p <- betas %>% ggplot(aes(x=Difference,y=Effect,label=Filter)) +
      geom_point(aes(size=Information),color=ifelse(abs(betas$Effect)>.01 & abs(betas$Difference)>.01,'#00BFC4','grey50')) +
      geom_text_repel(aes(label=ifelse(abs(Effect)>.01 & abs(Difference)>.01,as.character(Filter),'')),hjust=0,vjust=0,box.padding=.3)
  } else {
    p <- betas %>% ggplot(aes(x=Difference,y=Effect,label=Filter)) +
      geom_point(aes(size=Information),color=ifelse(abs(betas$Effect)>.01 & abs(betas$Difference)>.01 & grepl('de Novo',betas$Filter),'#F8766D',ifelse(abs(betas$Effect)>.01 & abs(betas$Difference)>.01 & !grepl('de Novo',betas$Filter),'#00BFC4','grey50'))) +
      geom_text_repel(aes(label=ifelse(abs(Effect)>.01 & abs(Difference)>.01,as.character(Filter),'')),hjust=0,vjust=0,box.padding=1)
  }
  p <- p +
    geom_hline(yintercept=0,col='black',lty=2) +
    geom_vline(xintercept=0,col='black',lty=2) +
    scale_size(guide = FALSE) +
    theme(panel.background = element_rect(fill = NA, color = "black"),
          axis.title=element_text(size=12),
          axis.text=element_text(size=12)) +
    ylab('Effect Size') +
    xlab('Mean Activation Difference')
  if (!is.null(fl)) {ggsave(fl,p,width=15,height=7.5)} else {return(p)}
}
#' Plot first-layer filter activations by class.
#'
#' Plot the first-layer convolutional filter activations by class label and filter.
#'
#' @param data Tibble containing effect size information, filter name, mean activation difference, and information.
#' @param show_beta Whether to include the estimated effect size estimates as well.
#' @param show_offset Whether to include the estimated filter offsets as well.
#' @param fl Filename to save plots to. If not provided, a plot object will be returned.
#'
#' @return A ggplot object.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords plot activations filter
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_boxplot geom_hline theme labs ylab guides ggsave element_rect element_text guide_legend
#' @export
plot_filter_activations_byClass <- function(data,show_beta=F,show_offset=F,fl=NULL) {
  if ('Location' %in% colnames(data)) {
    data$Location <- factor(data$Location, levels = order(unique(as.numeric(data$Location))))
    data$Filter <- as.factor(paste('Filter',data$Filter))
    if (show_beta | show_offset) {
      p <- ggplot(data,aes(x=Location,y=Activation,fill=Y),show.legend=TRUE) +
        geom_boxplot()
      if (show_beta) {
        p <- p + geom_point(aes(x=Location,y=W_fc1),color='red',shape=4, size=2.5,show.legend = F)
      }
      if (show_offset) {
        p <- p + geom_point(aes(x=Location,y=b_conv1),color='blue',shape=3, size=2.5,show.legend = F)
      }
    } else {
      p <- ggplot(data,aes(x=Location,y=Activation,fill=Y)) +
        geom_boxplot
    }
    p <- p +
      facet_wrap(vars(Filter)) +
      geom_hline(yintercept=0,col='black') +
      geom_hline(yintercept=1,col='black') +
      geom_hline(yintercept=.5,lty=4,col='black') +
      labs(title='First-layer Filter Activations by Location') +
      theme(panel.background = element_rect(fill = NA, color = "black"),
            axis.title=element_text(size=20),
            axis.text=element_text(size=14),
            title = element_text(size=24),
            legend.title = element_text(size=14),
            legend.text = element_text(size=14),
            legend.position = 'right',
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab('P(Z=1)') +
      geom_point(aes(shape='Beta_1'),alpha=0) +
      geom_point(aes(color='Beta_0'),alpha=0)
    if (show_beta & show_offset) {
      p <- p + guides(fill=guide_legend(order=1),
                      shape=guide_legend(title=NULL, override.aes = list(alpha = 1,color='red',shape=4,size=2.5,order=2)),
                      color=guide_legend(title=NULL, override.aes = list(alpha = 1,color='blue',shape=3,size=2.5,order=3)))
    } else if (show_beta) {
      p <- p + guides(fill=guide_legend(order=1),
                      shape=guide_legend(title=NULL, override.aes = list(alpha = 1,color='red',shape=4,size=2.5,order=2)))
    } else if (show_offset) {
      p <- p + guides(fill=guide_legend(order=1),
                      color=guide_legend(title=NULL, override.aes = list(alpha = 1,color='blue',shape=3,size=2.5,order=2)))
    } else {
      p <- p + guides(fill=guide_legend(order=1))
    }
  } else {
    if (show_beta | show_offset) {
      p <- ggplot(data,aes(x=Filter,y=Activation,fill=Y),show.legend=TRUE) +
        geom_boxplot()
      if (show_beta) {
        p <- p + geom_point(aes(x=Filter,y=W_fc1),color='red',shape=4, size=2.5,show.legend = F)
      }
      if (show_offset) {
        p <- p + geom_point(aes(x=Filter,y=b_conv1),color='blue',shape=3, size=2.5,show.legend = F)
      }
    } else {
      p <- ggplot(data,aes(x=Filter,y=Activation,fill=Y)) +
        geom_boxplot
    }
    p <- p +
      geom_hline(yintercept=0,col='black') +
      geom_hline(yintercept=1,col='black') +
      geom_hline(yintercept=.5,lty=4,col='black') +
      labs(title='First-layer Filter Activations') +
      theme(panel.background = element_rect(fill = NA, color = "black"),
            axis.title=element_text(size=20),
            axis.text=element_text(size=14),
            title = element_text(size=24),
            legend.title = element_text(size=14),
            legend.text = element_text(size=14),
            legend.position = 'right',
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylab('P(Z=1)') +
      geom_point(aes(shape='Beta_1'),alpha=0) +
      geom_point(aes(color='Beta_0'),alpha=0)
    if (show_beta & show_offset) {
      p <- p + guides(fill=guide_legend(order=1),
                      shape=guide_legend(title=NULL, override.aes = list(alpha = 1,color='red',shape=4,size=2.5,order=2)),
                      color=guide_legend(title=NULL, override.aes = list(alpha = 1,color='blue',shape=3,size=2.5,order=3)))
    } else if (show_beta) {
      p <- p + guides(fill=guide_legend(order=1),
                      shape=guide_legend(title=NULL, override.aes = list(alpha = 1,color='red',shape=4,size=2.5,order=2)))
    } else if (show_offset) {
      p <- p + guides(fill=guide_legend(order=1),
                      color=guide_legend(title=NULL, override.aes = list(alpha = 1,color='blue',shape=3,size=2.5,order=2)))
    } else {
      p <- p + guides(fill=guide_legend(order=1))
    }
  }
  if (!is.null(fl)) {ggsave(fl,p,width=15,height=7.5)} else {return(p)}
}
#' Plot filter motifs.
#'
#' Plot first-layer convolutional filter motifs.
#'
#' @param W_conv Either a list of motifs as ICMs, or a single ICM to visualize as sequence logos.
#' @param ylow Minimum Y-axis limit.
#' @param yhigh Maximum Y-axis limit.
#' @param method Y-label calculation type for sequence logo plot. Either \code{'custom'} or \code{'bits'}.
#' @param plotheight Height in inches for output plot.
#' @param fl Filename to save plots to. If not provided, a plot object will be returned.
#'
#' @return A ggplot object.
#'
#' @author Matthew Ploenzke, \email{ploenzke@@g.harvard.edu}
#' @keywords plot motif filter
#'
#' @importFrom ggplot2 ggplot aes geom_hline facet_wrap ylim labs theme ggsave element_rect element_text
#' @importFrom ggseqlogo geom_logo theme_logo
#' @export
plot_motifs <- function(W_conv,ylow=0,yhigh=2,method='custom',plotheight=NULL,fl) {
  if (is.null(plotheight)) {plotheight <- ifelse(opt$n_filters>8,15,7.5)}
  if (is.list(W_conv)) {
    p <- ggplot()+
      geom_logo(lapply(W_conv,function(i) {if (all(i==0)) {i+1e-10} else {i}}),method=method,seq_type='dna') +
      theme_logo() +
      geom_hline(yintercept=0,size=.25) +
      facet_wrap(~seq_group, ncol=4, scales= 'free_x')
  } else {
    p <- ggplot()+
      geom_logo(W_conv,method=method,seq_type='dna') +
      theme_logo() +
      geom_hline(yintercept=0,size=.25) +
      ylim(low=ylow,high=yhigh)
  }
  p <- p +
    labs(y='Bits') +
    theme(panel.background = element_rect(fill = NA, color = "black"),
          strip.text.x = element_text(size = 20),
          axis.title=element_text(size=20),
          axis.text=element_text(size=14),
          axis.text.x=element_text(size=14))
  if (!is.null(fl)) {ggsave(fl,p,width=15,height=plotheight,units='in')} else {return(p)}
}

