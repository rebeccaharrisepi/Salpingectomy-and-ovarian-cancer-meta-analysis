
# Packages used
library(readxl)
library(tidyverse)
library(metafor)

########################
### Read and set up data
########################
dat <- read_excel("Data.xlsx")

dat <- dat %>%
  filter(!Author=="Harnod", !Author== "Hanley") %>%
  mutate(logOR=log(Estimate),
         logUL=log(CIUL),
         logLL=log(CILL),
         SE= (logUL-logLL)/3.92,
         ESID = row_number(),
         sample = 1/Total) %>%
  select(Study, ESID, Author, Year, Location, Total, 'Study Design', Age_group, Time_sincesurgery, OC_type, salping_type, Effect, logOR, SE, Intervention_subgroups)

dat <- escalc(yi=logOR, sei=SE, data=dat)

order <- c("Cohort", "Nested Case-control", "Case-control")

dat <- dat %>%
  mutate(Study.Design = factor(Study.Design, 
                               levels = order,
                               labels = c("Cohort", "Nested Case-control","Case-control"))) %>%
  arrange(Study.Design, Year)


###########################################################
### Overall analysis (bilateral or unilateral salpingectomy)
###########################################################

data <- dat %>%
  filter(Intervention_subgroups==1)

overall <- rma(yi, vi=vi, subset=salping_type=="Bilateral or Unilateral", test="hksj",
               method="REML", slab=paste(Author, Year), ni=Total, data=data)
weights <- weights(overall)
mlabfun <- function(text, x) {
  list(bquote(paste(bold(.(text)),
                    "     (",tau^2, " = ", .(formatC(x$tau2, digits=3, format="f")), ", ",
                    I^2, " = ", .(formatC(x$I2, digits=1, format="f")), "%",
                    ")")))}

# Figure 2
forest(overall, atransf=exp, ilab=cbind(as.character(data$Location), as.character(data$Study.Design)),
       ilab.lab=c("Country", "Study Design"), ilab.xpos=c(-3.5, -2.5),
       predstyle="bar", 
       mlab=mlabfun("Bilateral or Unilateral Salpingectomy", overall))



####################################################
### Subgroup analysis based on type of salpingectomy
####################################################
data <- data %>%
  arrange(salping_type, Study.Design, Year)

Bilateral <- rma(yi, vi=vi, subset=salping_type=="Bilateral", 
                 method="REML", slab=paste(Author, Year), ni=Total, data=data)
Unilateral <- rma(yi, vi=vi, subset=salping_type=="Unilateral", 
                  method="REML", slab=paste(Author, Year), ni=Total, data=data)

###Figure 3
data <- data %>%
  filter(salping_type=="Bilateral"| salping_type=="Unilateral")
forest(data$yi, data$vi, rows=c(21:15, 8:4), ylim=c(0,26), atransf=exp,
       slab=paste(data$Author, data$Year), 
       ilab=cbind(data$Location, as.character(data$Study.Design)), ilab.lab=c("Country", "Study Design"), ilab.xpos=c(-8.5, -6))
text(-13.3, c(22.5,9.5), pos=4, c("Bilateral Salpingectomy",
                                   "Unilateral Salpingectomy"), font=4)
addpoly(Bilateral, mlab=mlabfun("Pooled effect", Bilateral), predstyle="bar", row=13.5)
addpoly(Unilateral, mlab=mlabfun("Pooled effect", Unilateral), predstyle="bar", row=2.5)


# multilevel meta-analysis to compare unilateral and bilateral subgroups 
rho <- 0.6 ### correlation for between outcomes
V1 <- vcalc(vi=vi, cluster=Author, type=salping_type, rho=rho, data=data)
chemodel <- rma.mv(yi=yi, V=V1, mods=yi ~ 0+ salping_type, slab=paste(Author, Year), data=data, random = ~ 1 | Author/ESID, test="t", method = "REML")
chemodel


####################################################
### Subgroup analysis based on time since surgery
####################################################
###Time 
time <- dat %>%
  filter(Time_sincesurgery%in% c("<10", ">10") & salping_type=="Bilateral or Unilateral") %>%
  arrange(Time_sincesurgery)


lessthan10 <- rma(yi, vi, data=time, subset=Time_sincesurgery=="<10",
           method="REML")
lessthan10

morethan10 <- rma(yi, vi, data=time, subset=Time_sincesurgery==">10",
                  method="REML")
morethan10

# Figure 4
forest(time$yi, time$vi, rows=c(11:10, 5:3), ylim=c(0,15), atransf=exp,
       slab=paste(time$Author, time$Year), 
       ilab=cbind(time$Location, as.character(time$Study.Design)), ilab.lab=c("Country", "Study Design"), ilab.xpos=c(-2, -1.5))
addpoly(lessthan10, mlab=mlabfun("Pooled effect", lessthan10), predstyle="bar", row=9)
addpoly(morethan10, mlab=mlabfun("Pooled effect", morethan10), predstyle="bar", row=2)
text(-2.95, c(12,6), pos=4, c("<10 years since surgery",
                                  ">10 years since surgery"), font=4)


####################################################
### Subgroup analysis based on age at surgery
####################################################
###Age <50

res <- rma(yi, vi, data=dat, subset=Age_group=="<50",
           method="REML")
res
forest(res, atransf=exp, slab=paste(Author, Year),ilab=cbind(Location, as.character(Study.Design)), ilab.lab=c("Country", "Study Design"), ilab.xpos=c(-1.6, -1.2), mlab=mlabfun("Age <50 years", res), predstyle="bar")


####################
### Publication bias
####################

###Funnel plot
par(mar=c(5,6,4,1)+.1)
funnel(overall, yaxis="sei", level=c(90, 95, 99), shade=c("white", "gray55", "gray75"), back="gray90", hlines="gray90",
       digits=2, refline=0, atransf=exp, label=0, legend=T, ylab="Standard Error")

###Peter's test for funnel plot asymmetry
regtest(overall, model="lm", predictor="sei")









########################################################################
### process used to pool unilateral and bilateral results within studies 
### to obtain overall effect sizes
########################################################################


Dareliusuni <- dat %>%
  filter(Author=="Darelius" & salping_type=="Unilateral" & OC_type %in% c("Type 1", "Type 2"))
res <- rma(yi, vi, method="FE", data=Dareliusuni)
forest(res, atransf=exp, slab=Dareliusuni$OC_type)


Dareliusbi <- dat %>%
  filter(Author=="Darelius" & salping_type=="Bilateral" & OC_type %in% c("Type 1", "Type 2"))
res <- rma(yi, vi, method="FE", data=Dareliusbi)
forest(res, atransf=exp, slab=Dareliusbi$OC_type)

Dareliusunibioruni <- dat %>%
  filter(Author=="Darelius" & salping_type=="Bilateral or Unilateral" & OC_type %in% c("Type 1", "Type 2"))
res <- rma(yi, vi, method="FE", data=Dareliusunibioruni)
forest(res, atransf=exp, slab=Dareliusunibioruni$OC_type)


Madsen <- dat %>%
  filter(Author=="Madsen" & salping_type %in% c("Unilateral", "Bilateral"))
res <- rma(yi, vi, method="FE", data=Madsen)
forest(res, atransf=exp, slab=Madsen$salping_type)


Duus <- dat %>%
  filter(Author=="Duus" & salping_type %in% c("Unilateral", "Bilateral") & OC_type=="Epithelial")
res <- rma(yi, vi, method="FE", data=Duus)
forest(res, atransf=exp, slab=salping_type)


### Duus combining 10-19 and 20+ years
Duus <- dat %>%
  filter(Author=="Duus" & Time_sincesurgery %in% c("10–19", ">20"))
res <- rma(yi, vi, method="FE", data=Duus)
forest(res, atransf=exp, slab=Time_sincesurgery)


Falconer <- dat %>%
  filter(Author=="Falconer" & Time_sincesurgery %in% c("0–4", "5–9") & salping_type=="Bilateral or Unilateral")
res <- rma(yi, vi, method="FE", data=Falconer)
forest(res, atransf=exp, slab=Time_sincesurgery)



