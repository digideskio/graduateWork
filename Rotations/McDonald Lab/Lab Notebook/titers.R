
# Original RIO3 and NSV experiments
# Probs not what we want to ultimately use...
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
        sd(data[data$kinase=="RIO3",]$PFU)/sqrt(5))

par(mar = c(5, 6, 4, 5))
plotTop <- max(1e+9)
barCenters <- barplot(height = mean, names.arg=names,
                      col=c("darkslategray", "gray80"), space=0.5,
                      las=1, ylim=c(1e+06,plotTop), cex.names = 0.75,
                      main = "Rotaviral titers following LV transduction",
                      ylab = "", xlab = "", border = "black",
                      axes = TRUE, log = "y", yaxt="n")
y=c(1e+06,1e+07,1e+08,1e+09)
ylab=c("1x10^6","1x10^7","1x10^8","1x10^9")
axis(2, at=y, labels=ylab, cex.axis = 0.75, las = 1)
mtext("PFU / mL", side=2,
      line=4.5, cex.lab=2, las=0, col="black")
mtext("Lentiviral vector", side=1,
      line=3, cex.lab=2, las=1, col="black")
segments(barCenters, mean-se*2, barCenters, mean+se*2, lwd=1.5)
arrows(barCenters, mean-se*2, barCenters, mean+se*2, lwd=1.5,
       angle=90, code=3, length = 0.05)

# install.packages("devtools")
library("devtools")
# install_github("ropensci/plotly")
library(plotly)
set_credentials_file("faulconbridge", "5vvex1bdyh")
# Fill in with your personal username and API key
# or, use this public demo account
py <- plotly()

trace1 <- list(
  x = c("NSV", "RIO3"), 
  y = c(1.9e+08, 2.6e+08), 
  name = "LVs", 
  error_y = list(
    type = "data", 
    array = c(24494897, 79686887),
    visible = TRUE
  ), 
  type = "bar",
  color = "rgba(31, 119, 180, 0.75)"
)

data <- list(trace1)
layout <- list(barmode = "group",
               autosize = FALSE, 
               width = 900, 
               height = 550, 
               margin = list(
                 l = 75, 
                 r = 25, 
                 b = 75, 
                 t = 50, 
                 pad = 5
               ), 
               paper_bgcolor = "#ffffff",
               plot_bgcolor = "#fcfcfc",
               title = "Rotaviral titers following LV transduction",
               titlefont = list(
                 famile = "Arial, sans-serif",
                 size = 30,
                 color = "black"
               ),
               xaxis = list(
                 type = "category", 
                 autorange = TRUE,
                 title = "Lentiviral vector", 
                 titlefont = list(
                   family = "Arial, sans-serif", 
                   size = 24, 
                   color = "black"
                 ), 
                 showticklabels = TRUE, 
                 tickangle = 0, 
                 tickfont = list(
                   family = "Arial, sans-serif", 
                   size = 14, 
                   color = "black"
                 ), 
                 exponentformat = "e", 
                 showexponent = "None",
                 showline = "true"
               ), 
               yaxis = list(
                 type = "log", 
                 autorange = TRUE,
                 title = "PFU / mL",
                 titlefont = list(
                   family = "Arial, sans-serif",
                   size = 24,
                   color = "black"
                 ),
                 showticklabels = TRUE,
                 tickangle = 0,
                 tickfont = list(
                   family = "Arial, sans-serif",
                   size = 14,
                   color = "black"
                 ),
                 showline = "true",
                 exponentformat = "power",
                 showexponent = "All"
               ))
response <- py$plotly(data, kwargs=list(layout=layout,
                                        filename="screen1",
                                        fileopt="overwrite"))
url <- response$url

########################################
# Titering: round 2!!!
########################################

data2 <- data.frame(kinase = c(rep("NSV",3),rep("RIO3",3),
                               rep("TGFa",3),rep("SKAP1",3),
                               rep("EPHA1",3),rep("CSNK2B",3),
                               rep("GUCY2D",3)),
                   vial = c(rep(seq(from = 1, to = 3, by = 1),7)),
                   PFU = c(400000000, 100000000, 200000000,
                           100000000, 100000000, 150000000,
                           300000000, 175000000, 300000000,
                           100000000, 100000000, 200000000,
                           250000000, 500000000, 200000000,
                           100000000, 200000000,  80000000,
                           100000000, 100000000, 400000000))
data2$kinase <- factor(data2$kinase)

names <- c("NSV", "RIO3", "TGFa", "SKAP1", "EPHA1", "CSNK2B", "GUCY2D")
mean <- c(mean(data2[data2$kinase=="NSV",]$PFU),
          mean(data2[data2$kinase=="RIO3",]$PFU),
          mean(data2[data2$kinase=="TGFa",]$PFU),
          mean(data2[data2$kinase=="SKAP1",]$PFU),
          mean(data2[data2$kinase=="EPHA1",]$PFU),
          mean(data2[data2$kinase=="CSNK2B",]$PFU),
          mean(data2[data2$kinase=="GUCY2D",]$PFU))
se <- c(sd(data2[data2$kinase=="NSV",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="RIO3",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="TGFa",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="SKAP1",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="EPHA1",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="CSNK2B",]$PFU)/sqrt(3),
        sd(data2[data2$kinase=="GUCY2D",]$PFU)/sqrt(3))

par(mar = c(5, 6, 4, 1))
plotTop <- max(1e+9)
barCenters <- barplot(height = mean, names.arg=names,
                      col=c("darkslategray", "gray80", "black",
                            "gray40", "gray60", "white", "darkgray"),
                      space=0.25,
                      las=1, ylim=c(1e+06,plotTop), cex.names = 0.7,
                      main = "Rotaviral titers following LV transduction",
                      ylab = "", xlab = "", border = "black",
                      axes = TRUE, log = "y", yaxt="n")
y=c(1e+06,1e+07,1e+08,1e+09)
ylab=c("1x10^6","1x10^7","1x10^8","1x10^9")
axis(2, at=y, labels=ylab,
     cex.axis = 0.75, las = 1)
mtext("PFU / mL", side=2,
      line=4, cex.lab=2, las=0, col="black")
mtext("Lentiviral vector", side=1,
      line=3, cex.lab=2, las=1, col="black")
segments(barCenters, mean-se, barCenters, mean+se, lwd=1.5)
arrows(barCenters, mean-se, barCenters, mean+se, lwd=1.5,
       angle=90, code=3, length = 0.05)

trace1 <- list(
  x = c("NSV", "RIO3", "TGF-&alpha;", "SKAP1",
        "GUCY2D", "CSNK2B", "EPHA1"), 
  y = c(233333333, 116666667, 258333333, 1.02e+08,
        1.04e+07, 3.21e+08, 4.75e+08), 
  name = "LVs", 
  error_y = list(
    type = "data", 
    array = c(88191710, 16666667, 41666667, 29016524,
              4387512, 18265309, 47615234),
    visible = TRUE
  ), 
  type = "bar",
  color = "rgba(31, 119, 180, 0.75)"
)

data <- list(trace1)
layout <- list(barmode = "group",
               autosize = FALSE, 
               width = 900, 
               height = 550, 
               margin = list(
                 l = 75, 
                 r = 25, 
                 b = 75, 
                 t = 50, 
                 pad = 5
               ), 
               paper_bgcolor = "#ffffff",
               plot_bgcolor = "#fcfcfc",
               title = "Rotaviral titers following LV transduction",
               titlefont = list(
                 famile = "Arial, sans-serif",
                 size = 30,
                 color = "black"
               ),
               xaxis = list(
                 type = "category", 
                 autorange = TRUE,
                 title = "Lentiviral vector", 
                 titlefont = list(
                   family = "Arial, sans-serif", 
                   size = 24, 
                   color = "black"
                 ), 
                 showticklabels = TRUE, 
                 tickangle = 0, 
                 tickfont = list(
                   family = "Arial, sans-serif", 
                   size = 14, 
                   color = "black"
                 ), 
                 exponentformat = "e", 
                 showexponent = "None",
                 showline = "true"
               ), 
               yaxis = list(
                 type = "log", 
                 autorange = TRUE,
                 title = "PFU / mL",
                 titlefont = list(
                   family = "Arial, sans-serif",
                   size = 24,
                   color = "black"
                  ),
                 showticklabels = TRUE,
                 tickangle = 0,
                 tickfont = list(
                   family = "Arial, sans-serif",
                   size = 14,
                   color = "black"
                 ),
                 showline = "true",
                 exponentformat = "power",
                 showexponent = "All"
               ))
response <- py$plotly(data, kwargs=list(layout=layout,
                                        filename="screen2",
                                        fileopt="overwrite"))
url <- response$url
