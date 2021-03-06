library(ggplot2)
library(grid)
library(RCurl)
library(ggpubr)

ff.dir <- getURL("https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/chrom_states_descr.txt")
chrom_states <- read.delim(text=ff.dir, as.is = TRUE)

element_custom <- function(...) {
  structure(list(...), class = c("element_custom", "element_blank"))
}

element_grob.element_custom <- function(element, label, x, y, ...)  {
  tg <- textGrob(label, y=y, gp=gpar(col=element$colour))
  padding <- unit(1,"line")
  rg <- rectGrob(y=y,width=grobWidth(tg)+padding, 
                 height=unit(1,"line")+padding, 
                 gp=gpar(fill = element$fill, col=NA, alpha=1))
  gTree(children=gList(rg, tg), width=grobWidth(tg) + padding, cl="custom_axis")
}

widthDetails.custom_axis <- function(x) x$width + unit(1,"mm") # fudge





plotChromStates <- function(x, chrom_states, state.x="state", 
                            tile=TRUE, tit, max.p=4){
  df <- merge(chrom_states, x, by.x="State", by.y=state.x)
  levels(df$beta) <- c("OR<1", "OR>1")
  df$values <- NA
  df$values[df$beta=="OR<1" & df$p.value<0.05/15] <- 1
  df$values[df$beta=="OR>1" & df$p.value<0.05/15] <- 2
  df$values <- as.factor(df$values)
  
  df$state2 <- paste0(df$Description, " (",
                      df$State, ")")
  chrom_states$state2 <- paste0(chrom_states$Description, " (",
                                chrom_states$State, ")")
  ll <- rev(chrom_states$state2[chrom_states$Order])
  df$State <- factor(df$state2, levels=ll)
  
  col2hex <- function(col, alpha) rgb(t(col2rgb(col)), 
                                      alpha=alpha, 
                                      maxColorValue=255)
  
  mycol <- rev(col2hex(chrom_states$color)[chrom_states$Order])
  
  if (!tile){
    df$log10pval <- -log10(df$p.value)
    df$log10pval[df$log10pval>max.p] <- max.p
    ggplot(df, aes(lab.exposure, State, size = log10pval)) +
    geom_point(aes(col=values)) + 
    ylab("Blood chromatin states (ROADMAP)") +
    xlab("Exposures") + 
    scale_colour_manual(values=c("red", "blue"),
                        labels = c("OR<1", "OR>1"),
                        limits = c(1,2),
                        na.translate = FALSE) +
    scale_size_continuous(breaks=c(2:max.p),
                          labels=c(2:max.p),
                          limits = c(2, max.p),
                          range = c(2,max.p)) +
    labs(colour = "Enrichment") + 
    theme(axis.text.y = element_custom(fill=mycol),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust=1, size=14,
                                     vjust = 0.5),
          axis.title.x=element_text(size=16,face="bold")) +
    guides(size = guide_legend(title = "-log10(p-adj)")) +
      theme(axis.ticks.length = unit(5, "pt")) + 
    ggtitle(tit)
  }
  else {
    ggplot(df, aes(lab.exposure, State)) +
    geom_raster(aes(fill = values), hjust=0) +
    ylab("Blood chromatin states (ROADMAP)") +
    xlab("Exposures") +
    scale_fill_manual(values = c('red','blue'),
                      labels = c("OR<1", "OR>1"),
                      na.translate = FALSE) +
    labs(fill = "Enrichment") + 
    theme(axis.text.y = element_custom(fill=mycol),
          axis.title.y = element_text(size = 16, face = "bold"),
          axis.text.x = element_text(angle = 90, hjust=1),
          axis.ticks = element_blank()) +
    ggtitle(tit)
  }
}








