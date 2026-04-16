# RVC Nassau Grouper

devtools::install_github('jeremiaheb/rvc')
library(rvc)
library(tidyverse)
library(cmdstanr)
library(bayesplot)
library(loo)
library(kableExtra)
theme_set(theme_bw())

years_keys<-c(1999:2012, 2014, 2016, 2018, 2022, 2024)
years_tort<-c(1999, 2000, 2004, 2006, 2008, 2010, 2012, 2014, 2016, 2018, 2021, 2023, 2024)

all_data<-list()
for(yr in unique(c(years_keys, years_tort))) {
  tryCatch({
    dat <- getRvcData(years = yr, regions = c("FLA KEYS", "DRY TORT"))
    dat$sample_data$PRIMARY_SAMPLE_UNIT <- as.character(dat$sample_data$PRIMARY_SAMPLE_UNIT)
    all_data[[as.character(yr)]] <- dat$sample_data
  }, error = function(e) message(paste("Skipping", yr, ":", e$message)))
}

all_samples<-bind_rows(all_data)
table(all_samples$YEAR,all_samples$REGION)

all_surveys<-all_samples%>%
  distinct(YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR, REGION, 
           STRAT, PROT, HABITAT_CD, DEPTH, LAT_DEGREES, LON_DEGREES)

nassau_pos<-all_samples %>%
  filter(SPECIES_CD == "EPI STRI") %>%
  group_by(YEAR, PRIMARY_SAMPLE_UNIT, STATION_NR) %>%
  summarize(count = sum(NUM), .groups="drop") %>%
  filter(count > 0)

nassau <- all_surveys %>%
  left_join(nassau_pos, by=c("YEAR","PRIMARY_SAMPLE_UNIT","STATION_NR")) %>%
  mutate(count = replace_na(count, 0),
         year = factor(YEAR),
         protected = factor(PROT, levels=c(0,1,2), 
                            labels=c("Open","Partial","No-Take")),
         region = factor(REGION),
         habitat = factor(HABITAT_CD))

dim(nassau)
sum(nassau$count > 0)
mean(nassau$count == 0)

nassau$present <- ifelse(nassau$count > 0, 1, 0)
nassau$depthS <- (nassau$DEPTH - mean(nassau$DEPTH, na.rm=TRUE)) / sd(nassau$DEPTH, na.rm=TRUE)

nassau <- nassau %>%
  mutate(
    relief = factor(case_when(
      grepl("HR", HABITAT_CD) ~ "High",
      grepl("MR", HABITAT_CD) ~ "Mid",
      grepl("LR", HABITAT_CD) ~ "Low",
      TRUE ~ "Unknown"), levels=c("High","Mid","Low")),
    reef_morph = factor(case_when(
      grepl("CONT", HABITAT_CD) ~ "Continuous",
      grepl("ISOL", HABITAT_CD) ~ "Isolated",
      grepl("SPGR", HABITAT_CD) ~ "SpurGroove",
      grepl("RUBB", HABITAT_CD) ~ "Rubble",
      TRUE ~ "Other")))

table(nassau$protected, nassau$region)
table(nassau$relief)
table(nassau$reef_morph)
table(nassau$present)
str(nassau)

write.csv(nassau, "nassau_cleaned.csv", row.names=FALSE)

#Exploration

ggplot(nassau, aes(x=count)) +
  geom_histogram(binwidth=0.5) +
  ggtitle("Nassau Grouper Count Distribution")

ggplot(filter(nassau, count>0), aes(x=count)) +
  geom_histogram(binwidth=0.5) +
  ggtitle("Counts When Present")

nassau %>%
  group_by(YEAR, protected) %>%
  summarize(prop = mean(present), .groups="drop") %>%
  ggplot(aes(x=YEAR, y=prop, color=protected)) +
  geom_point() + geom_line() +
  labs(y="Proportion Present") +
  ggtitle("Occurrence by Year and Protection")

nassau %>%
  group_by(YEAR, region) %>%
  summarize(prop = mean(present), .groups="drop") %>%
  ggplot(aes(x=YEAR, y=prop, color=region)) +
  geom_point() + geom_line() +
  ggtitle("Occurrence by Year and Region")

ggplot(nassau, aes(x=habitat, y=count)) +
  stat_summary() +
  ggtitle("Mean Count by Habitat")

ggplot(nassau, aes(x=relief, y=count)) +
  stat_summary() +
  ggtitle("Mean Count by Relief")

ggplot(nassau, aes(x=reef_morph, y=count)) +
  stat_summary() +
  ggtitle("Mean Count by Reef Type")

ggplot(filter(nassau, count>0), aes(x=DEPTH)) +
  geom_histogram(binwidth=2) +
  ggtitle("Depth when present")

nassau %>%
  group_by(year) %>%
  summarize(mean=mean(count), var=var(count)) %>%
  ggplot(aes(x=mean, y=var)) +
  geom_point() +
  geom_abline(slope=1, color="red") +
  ggtitle("Variance vs Mean")


# Prior Pred. Checks

set.seed(42)
nsim <- 200
matrixPPC<-model.matrix(present ~ year + protected + region + depthS + relief + reef_morph - 1, 
                        data=nassau)
Ncoef<-ncol(matrixPPC)
prior_b<-matrix(rnorm(nsim * Ncoef, mean=0, sd=3), nrow=nsim, ncol=Ncoef)
prior_logit<-prior_b %*% t(matrixPPC)
prior_p<-plogis(prior_logit)

prior_p_long <- data.frame(
  sim = rep(1:nsim, each=nrow(nassau)),
  p = as.vector(t(prior_p)))

ggplot(prior_p_long, aes(x=p)) +
  geom_histogram(bins=50, fill="steelblue", alpha=0.7) +
  labs(x="Prior Pred Prob of Occurrence", y="Count",
       title="Prior Pred. Check (Distribution of Predicted Probabilites)")

prior_prop<-apply(prior_p, 1, function(row) mean(rbinom(length(row), 1, row)))

ggplot(data.frame(prop=prior_prop), aes(x=prop)) +
  geom_histogram(bins=30, fill="steelblue", alpha=0.7) +
  geom_vline(xintercept=mean(nassau$present), color="red", linewidth=1, linetype="dashed") +
  labs(x="Sim. Proportion of Surveys with Present", y="Count",
       title="Prior Predictive Check: Simulated Occur. Rates",
       subtitle="Red line = observed proportion (1.7%)")

# Binomial Models

binMod <- cmdstan_model("binomialMatrix.stan")
#STAN CODE: 
#data {
#  int N;
#  int Ncoef;
#  array[N] int Y;
#  matrix[N, Ncoef] xMatrix;
#}
#parameters {
#  vector[Ncoef] b;
#}
#transformed parameters {
#  vector[N] logitmu;
#  logitmu = xMatrix * b;
#}
#model {
#  b ~ normal(0, 3);
#  Y ~ bernoulli_logit(logitmu);
#}
#generated quantities {
#  vector[N] Yrep;
#  vector[N] LL;
#  for(i in 1:N) {
#    Yrep[i] = bernoulli_logit_rng(logitmu[i]);
#    LL[i] = bernoulli_logit_lpmf(Y[i] | logitmu[i]);
#  }
#}



# Model 1 = Interaction
matrix1<-model.matrix(present ~ year * protected + region - 1, data=nassau)
binData1<-list(N=nrow(nassau), Ncoef=ncol(matrix1), Y=nassau$present, xMatrix=matrix1)
binFit1<-binMod$sample(data=binData1, chains=4, parallel_chains=4, refresh=0)
binFit1$diagnostic_summary()
binFit1$summary(variables="b") %>% select(-mad,-median) %>% kable(digits=2)

# Model 2 = add Depth
matrix2<-model.matrix(present ~ year + protected + region + depthS - 1, data=nassau)
binData2<-list(N=nrow(nassau), Ncoef=ncol(matrix2), Y=nassau$present, xMatrix=matrix2)
binFit2<-binMod$sample(data=binData2, chains=4, parallel_chains=4, refresh=0)
binFit2$diagnostic_summary()
binFit2$summary(variables="b") %>% select(-mad,-median) %>% kable(digits=2)

# Model 3 = add Relief
matrix3 <- model.matrix(present ~ year + protected + region + relief - 1, data=nassau)
binData3 <- list(N=nrow(nassau), Ncoef=ncol(matrix3), Y=nassau$present, xMatrix=matrix3)
binFit3 <- binMod$sample(data=binData3, chains=4, parallel_chains=4, refresh=0)
binFit3$diagnostic_summary()
binFit3$summary(variables="b") %>% select(-mad,-median) %>% kable(digits=2)

# Model 4 = add Depth and Relief
matrix4 <- model.matrix(present ~ year + protected + region + depthS + relief - 1, data=nassau)
binData4 <- list(N=nrow(nassau), Ncoef=ncol(matrix4), Y=nassau$present, xMatrix=matrix4)
binFit4 <- binMod$sample(data=binData4, chains=4, parallel_chains=4, refresh=0)
binFit4$diagnostic_summary()
binFit4$summary(variables="b") %>% select(-mad,-median) %>% kable(digits=2)

# Model 4b: add Reef Morphology
matrix4b <- model.matrix(present ~ year + protected + region + depthS + relief + reef_morph - 1, 
                         data=nassau)
binData4b <- list(N=nrow(nassau), Ncoef=ncol(matrix4b), Y=nassau$present, xMatrix=matrix4b)
binFit4b <- binMod$sample(data=binData4b, chains=4, parallel_chains=4, refresh=0)
binFit4b$diagnostic_summary()
binFit4b$summary(variables="b") %>% select(-mad,-median) %>% kable(digits=2)


# Model Comparison

LL1<-binFit1$draws(format="matrix", variables="LL")
LL2<-binFit2$draws(format="matrix", variables="LL")
LL3<-binFit3$draws(format="matrix", variables="LL")
LL4 <-binFit4$draws(format="matrix", variables="LL")
LL4b<- binFit4b$draws(format="matrix", variables="LL")

waicTable <- data.frame(
  Model=c("1: Year*Prot + Region",
          "2: Year + Prot + Region + Depth",
          "3: Year + Prot + Region + Relief",
          "4: Year + Prot + Region + Depth + Relief",
          "4b: + Reef Morphology"),
  WAIC=c(waic(LL1)$estimates["waic","Estimate"],
         waic(LL2)$estimates["waic","Estimate"],
         waic(LL3)$estimates["waic","Estimate"],
         waic(LL4)$estimates["waic","Estimate"],
         waic(LL4b)$estimates["waic","Estimate"]),
  LOOIC=c(loo(LL1)$estimates["looic","Estimate"],
          loo(LL2)$estimates["looic","Estimate"],
          loo(LL3)$estimates["looic","Estimate"],
          loo(LL4)$estimates["looic","Estimate"],
          loo(LL4b)$estimates["looic","Estimate"]))
waicTable %>% kable(digits=2)

# Select best fixed-effects model
bestFit <- binFit4b
bestLL <- LL4b

# Post. Pred. Check

Draws <- bestFit$draws(variables="Yrep", format="matrix")
ppc_dens_overlay(y=nassau$present, yrep=Draws[1:25,]) +
  ggtitle("Posterior Predictive Check — Best Model")
rm(Draws); gc()


# Convergence

draws_b <- bestFit$draws(format="array", variables="b")
colnames(matrix4b)
n_year <- length(levels(nassau$year))
non_year_idx <- (n_year+1):ncol(matrix4b)

mcmc_dens_overlay(draws_b, pars=paste0("b[", non_year_idx, "]")) +
  ggtitle("Chain Density Overlays Best Model (full)")

mcmc_trace(draws_b, pars=paste0("b[", non_year_idx, "]")) +
  ggtitle("Trace Plots best model (full)")

# Drop one WAIC

# Drop protection
matrix5 <- model.matrix(present ~ year + region + depthS + relief + reef_morph - 1, data=nassau)
binData5 <- list(N=nrow(nassau), Ncoef=ncol(matrix5), Y=nassau$present, xMatrix=matrix5)
binFit5 <- binMod$sample(data=binData5, chains=4, parallel_chains=4, refresh=0)
LL5 <- binFit5$draws(format="matrix", variables="LL")

# Drop region
matrix6 <- model.matrix(present ~ year + protected + depthS + relief + reef_morph - 1, data=nassau)
binData6 <- list(N=nrow(nassau), Ncoef=ncol(matrix6), Y=nassau$present, xMatrix=matrix6)
binFit6 <- binMod$sample(data=binData6, chains=4, parallel_chains=4, refresh=0)
LL6 <- binFit6$draws(format="matrix", variables="LL")

# Drop relief
matrix7 <- model.matrix(present ~ year + protected + region + depthS + reef_morph - 1, data=nassau)
binData7 <- list(N=nrow(nassau), Ncoef=ncol(matrix7), Y=nassau$present, xMatrix=matrix7)
binFit7 <- binMod$sample(data=binData7, chains=4, parallel_chains=4, refresh=0)
LL7 <- binFit7$draws(format="matrix", variables="LL")

# Drop depth
matrix8 <- model.matrix(present ~ year + protected + region + relief + reef_morph - 1, data=nassau)
binData8 <- list(N=nrow(nassau), Ncoef=ncol(matrix8), Y=nassau$present, xMatrix=matrix8)
binFit8 <- binMod$sample(data=binData8, chains=4, parallel_chains=4, refresh=0)
LL8 <- binFit8$draws(format="matrix", variables="LL")

#model 4 already drops reef morph (already fit)

waicTable2 <- data.frame(
  Model=c("Full (4b: Prot+Region+Depth+Relief+Morph)",
          "Drop Protection", "Drop Region",
          "Drop Relief", "Drop Depth", "Drop Reef Morphology"),
  WAIC=c(waic(LL4b)$estimates["waic","Estimate"],
         waic(LL5)$estimates["waic","Estimate"],
         waic(LL6)$estimates["waic","Estimate"],
         waic(LL7)$estimates["waic","Estimate"],
         waic(LL8)$estimates["waic","Estimate"],
         waic(LL4)$estimates["waic","Estimate"]))
waicTable2$deltaWAIC <- waicTable2$WAIC - min(waicTable2$WAIC)
waicTable2 %>% kable(digits=2)

