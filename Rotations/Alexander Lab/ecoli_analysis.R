########################################
#                                      #
# Pipeline for analysis of E. coli	   #
# multi-antibiotic resistance data	   #
#                                      #
########################################

########################################
# NOTE: Here we assume that your data
# take the general form:
# 
# +------------+--------+--------+----------+
# | Resistance | Sample |  Date  | Transect |
# +------------+--------+--------+----------+
# |     4      |    1   |  d/m/y |    15    |
# +------------+--------+--------+----------+
# |     1      |    2   |  d/m/y |    15    |
# +------------+--------+--------+----------+
# |     3      |    3   |  d/m/y |    15    |
# +------------+--------+--------+----------+
#
# I.e., that each observation has its
# own unique row. Multiple observations
# for purposes of analysis should not
# be contained on the same row.
#
# `Resistance' refers to the number of
# antibiotics to which a sample is 
# found to be resistant.
########################################


# Import our data from a CSV file
#data <- read.csv("~/path/to/file/data.csv", header = TRUE)

########################################
# Generate bogus data for illustration
########################################

  # Generate transects 1 --- 49
  # (Odd numbers only)
	transect <- NULL
	for(i in 1:49) {
	  if(i %% 2 != 0) {
	    transect <- c(transect, rep(i,6))
	  }
	}
  
  # Generate landuse types
	landuse <- NULL
	landuse[transect >= 31] <- "park"
	landuse[transect >= 19 & transect <= 29] <- "town"
	landuse[transect <= 17] <- "mixed"

  # Generate dates with observations every
  # other week
	tmp <- c(seq.Date(from = as.Date("2011-06-23"),
	                to   = as.Date("2011-10-27"),
	                by = "2 weeks"),
	         seq.Date(from = as.Date("2011-11-03"),
	                  to   = as.Date("2012-03-15"),
	                  by = "2 weeks"))

	dates <- NULL
	for(i in 1:20) {
	  dates <- c(dates, rep(as.Date(tmp[i], format = '%y-%m-%d'),6))
	}
	dates <- c(rep(dates[1:60],25), rep(dates[61:120],25))
	dates <- as.Date(dates, origin="1970-01-01")

  # Generate random rainfall data
	rainfall <- NULL
	tmp <- c(sample(1:100, 40, replace = TRUE))
	for(i in 1:40) {
	  rainfall <- c(rainfall, rep(tmp[i], 75))
	}

  # Stitch everything together into a data frame
	data <- data.frame(Resistance = c(sample(0:10, 3000, replace = TRUE)),
	                   Sample     = c(rep(seq(from = 1, to = 6, by = 1), 500)),
	                   Date       = dates,
	                   Season     = c(rep("dry", 1500), 
	                                   rep("wet", 1500)),
	                   Transect   = c(rep(transect, 20)),
	                   Landuse    = c(rep(landuse,  20)),
	                   Rainfall   = rainfall)

# View our data
View(data)

########################################
# Recode seasons based on sampling date
########################################
data$Season <- NA

# Convert our dates from factors to datetime objects
str(data$Date)
data$Date <- as.date(data$Date, format='%d/%m/%y')

# Convert to wet / dry seasons
# Recode based on given date ranges in DD / MM / YYYY
data$Season[data$Date >= "2011-07-13" & data$Date <= "2011-10-20"] <- "dry"
data$Season[data$Date >= "2011-11-03" & data$Date <= "2012-03-06"] <- "wet"
data$Season[data$Date >= "2012-04-11" & data$Date <= "2012-04-25"] <- "dry"


########################################
# Recode transects to land use
########################################

# Construct a new column
data$Landuse <- NA

# Recode our transects
data$Landuse[data$TR >= 31] <- "park"
data$Landuse[data$TR >= 19 & TR <= 29] <- "town"
data$Landuse[data$TR <= 17] <- "mixed"


########################################
# Convert DVs to factors
########################################

data <- within(data, {
  Landuse  <- factor(Landuse)
  Season   <- factor(Season)
  Sample   <- factor(Sample)
})

########################################
# Reshape data wide to long
########################################

install.packages("reshape")
library(reshape)

Data$Average <- NULL
head(Data)

# Reshape data from wide to long
reshapedData <- melt(Data, id=c("Date", "Sample", "TR", "Season", "Landuse"))

# Create MDR factors
reshapedData$MDR <- NA
reshapedData$MDR[reshapedData$value >= 3] <- "MDR"
reshapedData$MDR[reshapedData$value > 0 & reshapedData$value < 3]  <- "Resistant"
reshapedData$MDR[reshapedData$value == 0] <- "Susceptible"

# Convert transect and MDR to factors
reshapedData <- within(reshapedData, {
  TR    <- factor(TR)
  MDR   <- factor(MDR)
})

########################################
# Construct lmer given:
#
# No. ABR samples    - Outcome
# E. coli sample No. - Random effect
# Rainfall           - Random effect
# Transect No.       - Fixed effect
# Land use category  - Fixed effect
# Date sampled       - Fixed effect
# Season sampled     - Fixed effect
# 
########################################

install.packages("lme4")
library(lme4)

# We will start by predicting MDR by
# transect, landuse, date, and season
# (Transect/landuse and date/season may
# may be redundant; individual coefficients
# can be cut where justified)
#
# We also specify random slopes for Sample
MDR.model <- lmer(value ~ Landuse + Season + 
                    (1 | Sample), reshapedData)

# Tukey's HSD pairwise test
library(multcomp)
glht(MDR.model, linfct = mcp(TR = "Tukey"))

# Or as an ANOVA:
MDR.aov <- aov(value ~ TR + season 
	           + Error(sample / (TR + season)), data = reshapedData)
