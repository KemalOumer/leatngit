## ------------------------------------------------------------------------
library(SUMMER)
library(ggplot2)
library(gridExtra)



## ------------------------------------------------------------------------
library(readstata13)
filename <- "D:/ETBR71DT/ETBR71FL.DTA"
births <- read.dta13(filename, generate.factors = TRUE)   

library(rgdal)
mapfilename <- "D:/map2016/shps/sdr_subnational_boundaries.shp"
geo <- readOGR(mapfilename, verbose = FALSE)
Amat <- getAmat(geo, geo$REGNAME)

#other way#
shape = readOGR(dsn="D:/map2016/shps", layer="sdr_subnational_boundaries")
Amat <- getAmat(shape, shape$REGNAME)


## ------------------------------------------------------------------------
loc <- readOGR("D:/ETGE71FL/ETGE71FL.shp", verbose = FALSE)
loc.dat <- data.frame(cluster = loc$DHSCLUST, long = loc$LONGNUM, lat = loc$LATNUM)
gps <- mapPoints(loc.dat, geo = geo, long = "long", lat = "lat", names = c("REGNAME"))
colnames(gps)[4] <- "region"
sum(is.na(gps$region))


## ------------------------------------------------------------------------
unknown_cluster <- gps$cluster[which(is.na(gps$region))]
gps <- gps[gps$cluster %in% unknown_cluster == FALSE, ]
births <- births[births$v001 %in% unknown_cluster == FALSE, ]
births <- merge(births, gps[, c("cluster", "region")], by.x = "v001", by.y = "cluster", all.x = TRUE)
births$v024 <- births$region


## ------------------------------------------------------------------------
dat <- getBirths(data = births, strata = c("v023"),  year.cut = seq(1990, 2020, by = 1)) 
dat <- dat[,c("v001", "v002",  "v024","time","age","v005","strata","died")]
colnames(dat) <- c("clustid", "id", "region","time","age","weights","strata","died")
years <- levels(dat$time)
head(dat)


## -----------------------------------------------------------
dat_infant = getBirths(data = births, month.cut = c(1, 12), strata=c("v023"))


## -----------------------------------------------------------
direct0 <- getDirect(births = dat, years = years, regionVar = "region",
                     timeVar = "time", clusterVar = "~clustid + id", ageVar = "age", weightsVar = "weights", geo.recode = NULL)


## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
fit0 <- fitINLA(data = direct0, geo = NULL, Amat = NULL, year_label = years, year_range = c(1990, 2020), rw = 2, m=1)


## ------------------------------------------------------------------------
out0 <- getSmoothed(fit0, year_range = c(1990, 2020), year_label = years)


## ------------------------------------------------------------------------
random.time <- getDiag(fit0, field = "time", year_label = years)
plot(random.time, is.subnational=FALSE) + facet_grid(~label) + ggtitle("Compare temporal random effects: National Model") + ylab("Random Effects")


## ------------------------------------------------------------------------

## ------------------------------------------------------------------------

g1 <- plot(out0, year_label = years, year_med = as.numeric(years), is.subnational=FALSE) + ggtitle("National Smoothing Model") + ylim(c(0, 0.17))

## ------------------------------------------------------------------------
fit0.sub <- fitINLA(data = direct0, geo = geo, Amat = Amat, year_label = years,
                    year_range = c(1990, 2019), rw = 2, type.st = 4, m = 1)


## ------------------------------------------------------------------------
out0.sub <- getSmoothed(fit0.sub, Amat = Amat, year_range = c(1990, 2019), year_label = years)


## ------------------------------------------------------------------------
plot(out0.sub, is.subnational=FALSE,  data.add = direct0, option.add = list(point = "mean", by = "survey"), color.add = "orange", per1000=TRUE) + facet_wrap(~region, ncol = 7) + theme(legend.position = "none") + scale_color_manual(values = rep("gray20", 
                                                                                                                                                                                                                                                    47))



## ------------------------------------------------------------------------
mapPlot(data = subset(out0.sub, years %in% c(1990, 1995, 2000, 2005, 2010, 2015)), geo = geo, variables=c("years"), values = c("median"), by.data = "region", by.geo = "REGNAME", is.long = TRUE, border = "gray80", size = 0.2, ncol = 4, per1000 = TRUE, legend.label = "U5MR")


## ------------------------------------------------------------------------
breaks.hatch <- seq(0, 900, len = 10)
hatchPlot(data = subset(out0.sub, years %in% c(1990, 1995, 2000, 2005, 2010, 2015)), geo = geo, variables=c("years"), values = c("median"), by.data = "region", by.geo = "REGNAME", lower = "lower", upper = "upper", is.long=TRUE, ncol = 4, hatch = "gray50", per1000 = TRUE, breaks.CI = breaks.hatch) 


## ------------------------------------------------------------------------
hatchPlot(data = subset(direct0, years %in% c( 1990, 1995, 2000, 2005, 2010, 2015,2019)), geo = geo, variables=c("years"), values = c("mean"), by.data = "region", by.geo = "REGNAME", lower = "lower", upper = "upper", is.long=TRUE, ncol = 4, hatch = "gray50", per1000 = TRUE, breaks.CI = breaks.hatch) 


## ------------------------------------------------------------------------
random.time <- getDiag(fit0.sub, field = "time", year_label = years)
random.space <- getDiag(fit0.sub, field = "space", Amat = Amat)
random.spacetime <- getDiag(fit0.sub, field = "spacetime", year_label = years, Amat = Amat)


## ------------------------------------------------------------------------
plot(random.time, is.subnational=FALSE) + facet_wrap(~label) + ggtitle("Compare temporal random effects") + ylab("Random Effects")

## ------------------------------------------------------------------------
mapPlot(random.space, geo = geo, by.data = "region", by.geo = "REGNAME", variables = "label", values = c("median"), ncol = 2, is.long = TRUE) 


## ------------------------------------------------------------------------
plot(random.spacetime, is.subnational=FALSE) + facet_wrap(~region, ncol = 7) + ylab("Random Effects")+ ggtitle("Compare space-time interaction random effects")


####END######





## ------------------------------------------------------------------------
results0 <- out0.sub[, c("region", "years", "median", "upper", "lower")]
results2 <- direct0[, c("region", "years", "mean", "upper", "lower")]
results0$type <- "Smooth_Direct"
results2$type <- "Direct"
colnames(results0)[3]  <- colnames(results2)[3] <- "Est"
results <- rbind(results0, results2)
results$type <- factor(results$type, levels = c("Direct", "Smooth_Direct"))
results <- results[results$region != "All", ]


## ------------------------------------------------------------------------
pos <- position_dodge(width = 0.2)
ggplot(results, aes(x = years, y = Est, color = type, ymin = lower, ymax = upper)) + geom_point(position=pos, alpha = 0.8, size = 0.5) + geom_errorbar(position=pos, size = 0.5, alpha = 0.8) + facet_wrap(~region, ncol = 7) + theme_bw() + theme(legend.position = "bottom") + scale_colour_manual(values=c("steelblue", "orange", "red")) 



## ------------------------------------------------------------------------
g <- NULL
for(i in 1:7){
  t <- c( 1990, 1995, 2000, 2005, 2010, 2015,2019)[i]
  selectyear <- subset(results, type == "Smooth_Direct" & years == t)
  ordered <- selectyear[order(selectyear$Est), "region"]
  results$region <- factor(results$region, levels = ordered)
  pos <- position_dodge(width = 0.8)
  g[[i]] <- ggplot(subset(results, years == t), aes(y = logit(Est), x = region, color = type, ymin = logit(lower), ymax = logit(upper))) + geom_point(position=pos,  size = 0.5) + geom_errorbar(position=pos, size = 0.8, width = 0, alpha = 0.8) + facet_grid(~years) + theme_bw() + theme(legend.position = "bottom") + ylab("Logit(U5MR)") + coord_flip() + scale_colour_manual(values=c("steelblue", "orange", "red")) 
}
grid.arrange(grobs = g, ncol = 4)


## ------------------------------------------------------------------------
library(tidyr)
range <- range(c(0, results$Est), na.rm = TRUE)
results.wide <- spread(results[, c(1,2,3,6)], type, Est) 
results.wide$Direct[is.na(results.wide$Direct)] <- 0

g1 <- ggplot(results.wide, aes(x = Direct, y = Smooth_Direct, color = region)) + geom_point(alpha = 0.5) + geom_abline(intercept = 0, slope = 1, color = "red") + xlim(range)+ ylim(range) + theme(legend.position = "none")
grid.arrange(g1, ncol = 2)
