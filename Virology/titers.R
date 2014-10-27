
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
        sd(data[data$kinase=="RIO3",]$PFU/sqrt(5)))


plotTop <- max(mean+se*2)
barCenters <- barplot(height = mean, names.arg=names,
                      col="gray", las=1, ylim=c(0,plotTop))
segments(barCenters, mean-se*2, barCenters, mean+se*2, lwd=2)
arrows(barCenters, mean-se*2, barCenters, mean+se*2, lwd=2,
       angle=90, code=3)

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
                 size = 24,
                 color = "black"
               ),
               xaxis = list(
                 type = "category", 
                 autorange = TRUE,
                 title = "Lentiviral vector", 
                 titlefont = list(
                   family = "Arial, sans-serif", 
                   size = 18, 
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
                   size = 18,
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

trace1 <- list(
  x = c("NSV", "RIO3", "TGF-&alpha;", "SKAP1"), 
  y = c(1.9e+08, 2.6e+07, 5.72e+07, 1.02e+08), 
  name = "LVs", 
  error_y = list(
    type = "data", 
    array = c(24494897, 7968688, 18367259, 29016524),
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
                 size = 24,
                 color = "black"
               ),
               xaxis = list(
                 type = "category", 
                 autorange = TRUE,
                 title = "Lentiviral vector", 
                 titlefont = list(
                   family = "Arial, sans-serif", 
                   size = 18, 
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
                   size = 18,
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
