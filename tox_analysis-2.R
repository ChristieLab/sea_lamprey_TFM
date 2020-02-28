install.packages('ggplot2')
library(ggplot2)

setwd("E:/Purdue/pmarinus/TFM/Fig-v4/Fig2/")

### data and code for Lethal Concentration Measures
LC_data <- read.csv("Lethal Concentration Results.csv")

LC_data$Trial <- factor(LC_data$Trial,
                        levels = c("Tox 3", "Tox 2", "Tox 4"))


#plot without 99% mort
no_99 <- subset(LC_data, Percent < 99 )

no_99$Trial <- factor(no_99$Trial,
                      levels = c("Tox 3", "Tox 2", "Tox 4"))

no_99$Percent <- factor(no_99$Percent,
                        levels = c("25", "50"))

pdf("E:/Purdue/pmarinus/TFM/Fig-v4/Fig2/no_99.pdf")
ggplot(no_99, aes(x = Trial, y = LD, ymin=(LD-X), ymax = (LD + X), colour = Population)) +
  geom_point(size=4, position=position_dodge(width=1)) +
  geom_errorbar(size = 1, position=position_dodge(width=1)) +
  scale_color_manual(values=c("navy", "cadetblue2"),labels = c( "Michigan", "Massachusetts")) +
  facet_wrap(~Percent, labeller = labeller(Percent = c("25"="LC25", "50"="LC50"))) +
  theme_classic(base_size = 28) +
  theme(axis.text.x = element_blank(), axis.text.y = element_text(size=(20), colour = "black")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0,10), expand = c(0, 0), breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10))
  ggsave("no_99.pdf", width = 3, height = 4, units = "in")
dev.off()
