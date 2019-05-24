library(RCurl)
x <- getURL("https://raw.githubusercontent.com/isglobal-brge/brgeUtils/master/data/chrom_states_descr.txt")
chrom_states <- read.delim(text=x, as.is = TRUE)

df <- data.frame(state=rep(chrom_states$State,2),
                 pvalue=c(runif(15, 0, 1),runif(15, 0, 1)),
                 beta=rep(c(-1,-1,-1, -1, rep(1, 11)),2),
                 expos = rep(c("A", "B"), each=15))


df2 <- merge(df, chrom_states, by.x="state", by.y="State")

df2$State <- factor(df2$state, 
                    levels=rev(chrom_states$State[chrom_states$Order]))

library(ggplot2)
df2$state2 <- paste0(df2$Description, " (", df2$state, ")")
ggplot(df2, aes(expos, state2, col=beta)) +
  geom_point(aes(size = pvalue)) + 
  theme(axis.text.y = 
          element_text(colour = 
          rev(col2hex(chrom_states$color)[chrom_states$Order])))

df2$color
colors()
