data <- data.frame(kinase = c(rep("NSV",5),rep("RIO3",5)),
                   vial = c(rep(seq(from = 1, to = 5, by = 1),2)),
                   PFU = c(250000000,200000000,200000000,200000000,100000000,
                           150000000,100000000,550000000,200000000,300000000),
                   colorScheme = c(rep("Monochromatic",5),rep("Colored",5)))
data$kinase <- factor(data$kinase)

names <- c("NSV", "RIO3")
mean <- c(mean(data[data$kinase=="NSV",]$PFU),
          mean(data[data$kinase=="RIO3",]$PFU))
se <- c(sd(data[data$kinase=="NSV",]$PFU)/sqrt(5),
        sd(data[data$kinase=="RIO3",]$PFU/sqrt(5)))


plotTop <- max(mean+se*2)
barCenters <- barplot(height = mean, names.arg=names,
                      col="gray", las=1, ylim=c(0,plotTop))
segments(barCenters, mean-se*2, barCenters, mean+se*2, lwd=2)
arrows(barCenters, mean-se*2, barCenters, mean+se*2, lwd=2,
       angle=90, code=3)
