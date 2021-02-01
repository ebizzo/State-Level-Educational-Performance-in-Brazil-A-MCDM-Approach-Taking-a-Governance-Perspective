################################################################################
## State level Educational Perfomance
## Authors: Jorge Antunes
##          Peter Wanke
##          Eduardo Bizzo
##          Rodrigo Furst 
##          Franklin Mixon 
## Last Update: 2021-01-31
################################################################################


################################################################################
## Clean Memory
################################################################################
rm(list=ls())
gc()


################################################################################
## Libraries
################################################################################
## Install Packages - if necessary
## install.packages(c("readxl", "writexl", "topsis", "DEoptim", 
##                    "entropy", "combinat", "lpSolve", "lpSolveAPI", 
##                    "ggplot2", "reshape2", "ggrepel", "ggpubr", "AER", 
##                    "betareg", "simplexreg", "doParallel", 
##                    "lmtest", "flexmix", "psych"), 
##                  dependencies = TRUE)


library(readxl)
library(writexl)
library(topsis)
library(DEoptim)
library(entropy)
library(combinat)
library(lpSolve)
library(lpSolveAPI)
library(ggplot2)
library(reshape2)
library(ggrepel)
library(ggpubr)
library(AER)
library(betareg)
library(simplexreg)
library(doParallel)
library(lmtest)
library(flexmix)
library(psych)


################################################################################
## Read Dataset
################################################################################
## Path
path <- "D:/Dropbox/State-Level Educational Performance in Brazil/Bases/"
setwd(path)

## Read dataset
rank <- read_excel(
  path = "Base State-Level Educational Performance in Brazil.xlsx"
) 


################################################################################
## Data Manipulation
################################################################################
## Correct names
colnames(rank) <- c("State",                                                                           
                    "Year",                                                                            
                    "Municipal Spending on Education and Culture",
                    "State Spending on Education and Culture",                                         
                    "IDEB Early Years",                                                                
                    "IDEB Final Years",                                                                
                    "IDEB High School",                                                                
                    "Literacy Rate",                                                                   
                    "School Attendance Rate",                                                          
                    "School Attainment Rate",                                                        
                    "North",                                                                           
                    "Northeast",                                                                       
                    "South",                                                                           
                    "Southeast",                                                                       
                    "Central-West",                                                                    
                    "GDP",                                                                             
                    "HDI",                                                                             
                    "Infant Mortality",                                                                
                    "% Enrolled in Pub. Schools")
colnames(rank) <- gsub(pattern = " ", replacement = "_", x = colnames(rank))

## Create a Region variable
rank$Region <- factor(
  ifelse(rank$North == 1, "North",
         ifelse(rank$Northeast == 1, "Northeast",
                ifelse(rank$South == 1, "South",
                       ifelse(rank$Southeast == 1, "Southeast",
                              "Central-West"))))
)


## Matrix
head(rank)

## Decision Matrix
d   <- as.matrix(rank[,3:10])
imp <- c("-","-","+","+","+","+","+","+")


################################################################################
## Maximal Entropy Weights for TOPSIS Scores
################################################################################
## Function to evaluate weights for maximal entropy topsis scores
fn <- function(w){
  ## Evaluate topsis
  topsis_res <- topsis(decision = d, weights = w, impacts = imp)
  
  ## Evaluate entropy
  entropy_res <- entropy(topsis_res$score)
  
  ## return entropy
  return(-entropy_res)
}

## Optimization (Maximize Entropy)
Optimal_weights_entropy <- DEoptim(fn, 
                                   lower = rep(0, ncol(d)), 
                                   upper = rep(1, ncol(d)),
                                   control = list(itermax = 2500, 
                                                  trace = 50,
                                                  reltol = 1e-16, 
                                                  steptol = 100,
                                                  parallelType = 1,
                                                  packages = list("entropy",
                                                                  "topsis"),
                                                  parVar = list("d", "imp")))


## Optimal weights evaluation (entropy)
w_opt_entropy <- Optimal_weights_entropy$optim$bestmem
w_opt_entropy <- w_opt_entropy/sum(w_opt_entropy)
names(w_opt_entropy) <- colnames(d)


################################################################################
## NG Weights for TOPSIS Scores
################################################################################
## Create a d_aux variable
d_aux <- d

## Sum all Spendings for feasibility
d_aux[,1] <- log(rowSums(exp(d_aux[,1:2])))

## sum all IDEBs
d_aux[,3] <- rowSums(d_aux[,3:5])

## Remove variabels IDEB and Despesas unused
d_aux <- d_aux[,-c(2,4,5)]

## New names
colnames(d_aux)[1:2] <- c("Desp", "IDEB")

## Create permutations with variables
permvar <- permn(colnames(d_aux))

## List with permutations and weights
listaperm <- list(comb = permvar, wvar = as.list(rep(NA, length(permvar)))) 

## Matrix with weights for each DMU and permutation
mat_weights <- data.frame(matrix(NA, 
                                 nrow = nrow(d_aux), 
                                 ncol = ncol(d_aux)*length(permvar)))

## Loop in all permutations
for (zz in 1:length(listaperm[[2]])){
  
  ## order y in same order of permutation
  y <- data.frame(d_aux[ ,c(listaperm$comb[[zz]]) ])
  
  product_number  <- nrow(y)  ## Number of DMUs
  criteria_number <- ncol(y)  ## Number of variables
  
  ## Matrices with decisions values
  w    <- matrix(0, nrow = product_number, ncol = criteria_number)
  wres <- matrix(0, nrow = product_number, ncol = criteria_number)
  Ind  <- matrix(0, nrow = product_number, ncol = 1)
  
  ## Matrix with input parameters 
  yn <- data.frame(y)
  
  ## Normalization between 0-1
  yn <- apply(yn, 2, function(x){x/max(x)})

  ## Loop in DMU to solve LP
  for (i in 1:product_number) {
    w[i,] <- 0
    lprec <- make.lp(0, criteria_number)
    set.objfn(lprec, c(yn[i,])) ## Objective function
    
    ## First Restriction
    add.constraint(lprec, c(w[i,]+1),"=",1) ## Sum to one
    
    ## Set of restrictions 2
    for (j in 1:(criteria_number-1)) {
      w[i,]<-0
      w[i,j]<- 1
      w[i,j+1]<- -1
      add.constraint(lprec,c(w[i,]),">=",0)
    }
    
    ## Set of restrictions 3
    for (j in 1:criteria_number) {
      w[i,]<-0
      w[i,j]<-1
      add.constraint(lprec,c(w[i,]),">=",0)
    }
    
    ## Solve linear programming
    lp.control(lprec,sense='max')
    solve(lprec)
    
    ## add results in matrices
    wres[i,]<-get.variables(lprec)
    Ind[i]<-get.objective(lprec)
  }
  
  ## Results of permutation in mat_weights
  mat_weights[,((zz-1)*ncol(d_aux) + 1):(((zz)*ncol(d_aux)))] <- (wres)
  colnames(mat_weights)[((zz-1)*ncol(d_aux) + 1):(((zz)*ncol(d_aux)))] <- 
    colnames(y) 
  
  ## Print Evaluation
  print(paste(round((100 * zz / length(listaperm[[2]])), digits = 2),
              " %", sep = ""))
}	

## Optimal weights (NG)
w_NG = data.frame(matrix(NA, ncol = ncol(d_aux), nrow = nrow(d_aux)))
colnames(w_NG) <- colnames(d_aux)

for(i in 1:ncol(w_NG)){
  w_NG[,i] <- rowMeans(mat_weights[,colnames(mat_weights) == colnames(w_NG)[i]])
}

## Mean Weights NG
w_NG <- colMeans(w_NG)

## rep weights for agregatted variables
w_NG <- w_NG[c(rep(1,2),rep(2,3),3:5)]

## Normalize to sum 1.0
w_NG <- w_NG/sum(w_NG)

## Correct names
names(w_NG) <- colnames(d)


################################################################################
## Evaluate all TOPSIS
################################################################################
## Topsis (NG weights)
topsis_NG <- topsis(decision = d, weights = w_NG, impacts = imp)

## Topsis for optimal weights (entropy weights)
topsis_entropy <- topsis(decision = d, 
                         weights = w_opt_entropy, 
                         impacts = imp)

## Topsis (equal weights)
topsis_eq  <- topsis(decision = d, 
                     weights = rep(1/ncol(d), ncol(d)), 
                     impacts = imp)


################################################################################
## Results
################################################################################
## Correlation between scores for 3 weights
setwd("../Results/")
data_scores<- data.frame(Equal = topsis_eq$score,
                         NG = topsis_NG$score,
                         Entropy = topsis_entropy$score)
colnames(data_scores) <- paste(colnames(data_scores), "Weights", sep = " ")
write_xlsx(x = data.frame(cor(data_scores)), 
           path = "Correlation Scores TOPSIS.xlsx")


## density plots
data_density <- melt(data_scores)

pdf(file = "Density Plot Scores.pdf", width = 6, height = 5)
ggplot(data = data_density) + 
  geom_density(mapping = aes(x = value, color = variable)) + 
  ylab("Density") + xlab("TOPSIS Score") +
  xlim(0,1.0) +
  scale_color_discrete("") + 
  theme_bw()
dev.off()

## Segregation in North/Northeast and South/Southeast/Central-West
data_region_scores <- data.frame(Region = as.character(data.frame(rank)$Region),
                                 data_scores[,-1])
colnames(data_region_scores) <- gsub(pattern = "\\.", 
                                     replacement = " ", 
                                     x = colnames(data_region_scores))
data_region_scores$Region <- 
  factor(ifelse(data_region_scores$Region %in% c("North", "Northeast"), 
                "North and Northeast",
                "South, Southeast, Central-West"))
data_region_scores <- melt(data_region_scores)

pdf(file = "Scores Region Segregation.pdf", width = 8, height = 5)
ggplot(data = data_region_scores) + 
  facet_wrap("Region", ncol = 2) + 
  geom_density(aes(x = value, lty = variable)) + 
  scale_color_discrete("") + 
  ylab("Density") + xlab("TOPSIS Score") +
  xlim(0,1.0) +
  theme_bw() +
  theme(legend.title = element_blank())
dev.off()


## Line plots of TOPSIS
data_line <- data.frame(data.frame(rank)[,c(1:2)],
                        data_scores[,-1])
colnames(data_line) <- c("DMU", "Year", "NG Weights", "Entropy Weights")

data_line_mean <- do.call("rbind", 
                          lapply(split(data_line, data_line$Year), 
                                 function(x){
                                   return(colMeans(x[,-c(1:2)]))
                                 }))

data_line_mean <- data.frame(DMU = "Mean", 
                             Year = as.numeric(rownames(data_line_mean)),
                             data_line_mean)
colnames(data_line_mean) <- colnames(data_line)
data_line <- rbind(data_line, data_line_mean)
rm(data_line_mean)

data_line <- melt(data_line, id.vars = c(1:2))

pdf(file = "TOPSIS line plot.pdf", width = 8, height = 5)
ggplot() + 
  facet_wrap("variable", ncol = 2) + 
  geom_line(aes(x = Year, y = value, fill = DMU), 
            show.legend = FALSE, 
            data = data_line[data_line$DMU %in% 
                               c("Pernambuco", "Ceará", "Mean") == FALSE, ], 
            color = "gray70") + 
  geom_line(aes(x = Year, y = value, color = DMU), 
            data = data_line[data_line$DMU %in% 
                               c("Pernambuco", "Ceará", "Mean"), ], lwd = 1) +
  geom_point(aes(x = Year, y = value, color = DMU), 
            data = data_line[data_line$DMU %in% 
                               c("Pernambuco", "Ceará", "Mean"), ], size = 2) +
  scale_color_discrete("") + 
  ylab("TOPSIS Score") + xlab("Year") + 
  ylim(0,1.0) + 
  theme_bw()
dev.off()

## Barplot of weights
names(w_opt_entropy) <- names(w_NG) <- 
  c("Spending on Education\nand Culture (Municipality)",
    "Spending on Education\nand Culture (State)",
    "IDEB Early Years",
    "IDEB Final Years",
    "IDEB High School",
    "Literacy Rate",
    "School Attendance Rate",
    "School Attainment Rate")

df_weights <- cbind(melt(w_opt_entropy), melt(w_NG))
df_weights <- data.frame(Variable = rownames(df_weights), df_weights)
colnames(df_weights) <- c("Variable", "Entropy", "NG")

df_weights <- melt(df_weights, id.vars = 1)
colnames(df_weights) <- c("Variables", "Type", "Weight")


pdf(file = "Optimal Weights.pdf", width = 7, height = 5)
ggplot(data = df_weights) + 
  geom_col(mapping = aes(x = Variables, y = Weight, fill = Type), 
           position = "dodge") +
  xlab("") + ylab("Weight") +
  scale_fill_discrete("") + 
  coord_flip() +
  theme_bw()
dev.off()

## Average TOPSIS Score per Region and State (Only NG Weights)
data_average <- data.frame(data.frame(rank)[,c(1, 20, 2)], data_scores[,2])
colnames(data_average) <- c("DMU", "Region", "Year", "NG Weights")

data_average <- do.call(rbind,
                        lapply(split(data_average, data_average$DMU), 
                               function(x){
                                 x_mean <- x[1,,drop = FALSE]
                                 x_mean[1,4] <- mean(x[,4])
                                 return(x_mean[,-3])
                               }))

pdf(file = "Average TOPSIS by Region.pdf", width = 9, height = 7)
ggplot(data = data_average) + 
  geom_point(aes(x = Region, y = `NG Weights`)) +
  geom_text_repel(aes(x = Region, y = `NG Weights`, label = DMU), 
                  seed = 42, max.iter = 5000, nudge_x = 0.1) +
  ylab("Average TOPSIS Score") + 
  theme_bw()
dev.off()


################################################################################
## Tobit Regression - Entire models
################################################################################
## Create dataset for regression
data_reg <- data.frame(Performance = topsis_NG$score,
                       rank[,c(11:19)])

## Change names
colnames(data_reg) <- c("Performance", "North", "Northeast", "South",
                        "Southeast", "Central-West", "GDP", "HDI", 
                        "Infant Mortality", "Enrolled Pub. Schools")


## Central-West as reference
data_reg$`Central-West` <- NULL

## Model 1 (Only Region Model)
frl_1 <- Performance ~ North + Northeast + South + Southeast

## Model 2 (Region and some Contextual Model)
frl_2 <- Performance ~ North + Northeast + South + Southeast +
  GDP + HDI + `Infant Mortality`

## Model 3 (Complete Model)
frl_3 <- Performance ~.

## Tobit regression
tobit_res_1 <- tobit(formula = frl_1, data = data_reg, left = 0, right = 1)
tobit_res_2 <- tobit(formula = frl_2, data = data_reg, left = 0, right = 1)
tobit_res_3 <- tobit(formula = frl_3, data = data_reg, left = 0, right = 1)

## Summary
summary(tobit_res_1)
summary(tobit_res_2)
summary(tobit_res_3)

## Check for Homoscedastic Residuals for all regressions
## Residuals plot
df_residual <- data.frame(Obs = 1:nrow(data_reg),
                          Model1 = resid(tobit_res_1),
                          Model2 = resid(tobit_res_2),
                          Model3 = resid(tobit_res_3))
colnames(df_residual) <- c("Obs", "Model 1", "Model 2", "Model 3")
df_residual <- melt(df_residual, id.vars = "Obs")

residual_plot <- ggplot(data = df_residual) + 
  facet_wrap("variable", ncol = 3, scales = "free") +
  geom_point(aes(x = Obs, y = value)) + 
  geom_smooth(aes(x = Obs, y = value), method = "gam") + 
  ylab("Resid") + xlab("Observation") +  
  ggtitle("Resid Plot") + 
  theme_bw()

## QQplot
df_qq <- data.frame(Model1 = resid(tobit_res_1),
                    Model2 = resid(tobit_res_2),
                    Model3 = resid(tobit_res_3))
colnames(df_qq) <- c("Model 1", "Model 2", "Model 3")
df_qq <- melt(df_qq)

plot_qq <- 
  ggplot(mapping = aes(sample = value), data = df_qq) + 
  facet_wrap("variable", ncol = 3, scales = "free") + 
  stat_qq() + 
  stat_qq_line() +
  ylab("Sample") + xlab("Theoretical") + 
  ggtitle("QQ plot") + 
  theme_bw()

## Check KL divergence of residuals with a normal sample
KL_residuals <- function(reg_resid){
  ## Required libraries
  require(ggplot2)
  require(reshape2)
  
  ## KL matrix
  kl_mat <- density(reg_resid, n = 1024)
  kl_mat <- data.frame(Quantiles = kl_mat$x,
                       Resid = kl_mat$y,
                       Normal = dnorm(x = kl_mat$x, 
                                      mean = 0, 
                                      sd = sd(reg_resid)))
  
  ## density plot
  kl_gg <- melt(kl_mat, id.vars = "Quantiles")
  Density_Plot <- 
    ggplot(data = kl_gg) +
    geom_line(aes(x = Quantiles, y = value, color = variable)) +
    ylab("Density") + xlab("Quantiles") +
    scale_color_discrete("") + 
    theme_bw()
  
  ## return list
  res <- list()
  res$DensityPlot <- Density_Plot
  res$Divergence <- KLdiv(as.matrix(kl_mat[,-1]))[1,2]
  
  ## return
  return(res)
}

kl_div_1 <- KL_residuals(reg_resid = resid(tobit_res_1))
kl_div_2 <- KL_residuals(reg_resid = resid(tobit_res_2))
kl_div_3 <- KL_residuals(reg_resid = resid(tobit_res_3))


## Density plot of residuals
density_plot_1 <- kl_div_1$DensityPlot + 
  ggtitle(paste("Model 1 - KL Divergence: ", 
                round(kl_div_1$Divergence, digits = 2),
                sep = ""))
density_plot_2 <- kl_div_2$DensityPlot + 
  ggtitle(paste("Model 2 - KL Divergence: ", 
                round(kl_div_2$Divergence, digits = 2),
                sep = ""))
density_plot_3 <- kl_div_3$DensityPlot + 
  ggtitle(paste("Model 3 - KL Divergence: ", 
                round(kl_div_3$Divergence, digits = 2),
                sep = "")) 

pdf(file = "Heteroskedasticity Residuals Check.pdf", width = 10, height = 8)
ggarrange(plotlist = list(residual_plot,
                          plot_qq,
                          ggarrange(plotlist = list(density_plot_1,
                                                    density_plot_2,
                                                    density_plot_3), 
                                    ncol = 3, nrow = 1)), 
          nrow = 3)
dev.off()


################################################################################
## Tobit Regression - Region models
################################################################################
## Create dataset for both regions
data_reg_N_NE    <- data_reg[rowSums(data_reg[,c(2:3)]) == 1, ]
data_reg_S_SE_CW <- data_reg[rowSums(data_reg[,c(2:3)]) == 0, ]

## Remove Region Dummies
data_reg_N_NE    <- data_reg_N_NE[,-c(2:5)]
data_reg_S_SE_CW <- data_reg_S_SE_CW[,-c(2:5)]

## Model Region 
frl_region <- Performance ~.

## Tobit regression
tobit_N_NE_res    <- tobit(formula = frl_region, data = data_reg_N_NE, 
                           left = 0, right = 1)

tobit_S_SE_CW_res <- tobit(formula = frl_region, data = data_reg_S_SE_CW, 
                           left = 0, right = 1)

## Summary
summary(tobit_N_NE_res)
summary(tobit_S_SE_CW_res)

## Check for Homoscedastic Residuals for all regressions
plot_region_reg_check <- function(tobit_res, region){
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  
  ## Residuals plot
  df_residual <- data.frame(Obs    = 1:length(resid(tobit_res)),
                            Resid  = resid(tobit_res))
  colnames(df_residual) <- c("Obs", region)
  
  df_residual <- melt(df_residual, id.vars = "Obs")
  
  residual_plot <- ggplot(data = df_residual) + 
    facet_wrap("variable", ncol = 3, scales = "free") +
    geom_point(aes(x = Obs, y = value)) + 
    geom_smooth(aes(x = Obs, y = value), method = "gam") + 
    ylab("Resid") + xlab("Observation") +  
    ggtitle("Resid Plot") + 
    theme_bw()
  
  ## QQ plot
  df_qq <- data.frame(Region = resid(tobit_res))
  colnames(df_qq) <- region
  df_qq <- melt(df_qq)
  
  plot_qq <- 
    ggplot(mapping = aes(sample = value), data = df_qq) + 
    facet_wrap("variable", ncol = 3, scales = "free") + 
    stat_qq() + 
    stat_qq_line() +
    ylab("Sample") + xlab("Theoretical") + 
    ggtitle("QQ plot") + 
    theme_bw()
  
  ## Check KL divergence of residuals with a normal sample
  kl_mat <- density(resid(tobit_res), n = 1024)
  kl_mat <- data.frame(Quantiles = kl_mat$x,
                       Resid = kl_mat$y,
                       Normal = dnorm(x = kl_mat$x, 
                                      mean = 0, 
                                      sd = sd(resid(tobit_res))))
  
  ## density plot
  kl_gg <- melt(kl_mat, id.vars = "Quantiles")
  Density_Plot <- 
    ggplot(data = kl_gg) +
    geom_line(aes(x = Quantiles, y = value, color = variable)) +
    ylab("Density") + xlab("Quantiles") +
    scale_color_discrete("") + 
    ggtitle(paste(region, "\nKL Divergence: ", 
                  round(KLdiv(as.matrix(kl_mat[,-1]))[1,2], digits = 2),
                  sep = "")) +
    theme_bw()
  
  ## Plot
  heteroske_plot <- ggarrange(
    plotlist = list(residual_plot,
                    plot_qq,
                    Density_Plot),
    ncol = 1, nrow = 3)
  
  ## return
  return(heteroske_plot)
  
}


pdf(file = "Heteroskedasticity Residuals Check Region Regression.pdf", 
    width = 9, height = 7)
ggarrange(
  plotlist = list(
    plot_region_reg_check(tobit_res = tobit_N_NE_res, 
                          region = "North & Northeast"),
    plot_region_reg_check(tobit_res = tobit_S_SE_CW_res, 
                          region = "Central-West, South & Southeast")
  ), 
  ncol = 2, nrow =1)
dev.off()


################################################################################
## Robustness Checks for Tobit Regression
################################################################################
## Beta regression
beta_res_1 <- betareg(formula = frl_1, data = data_reg)
beta_res_2 <- betareg(formula = frl_2, data = data_reg)
beta_res_3 <- betareg(formula = frl_3, data = data_reg)


## Simplex regression
simplex_res_1 <- simplexreg(formula = frl_1, data = data_reg)
simplex_res_2 <- simplexreg(formula = frl_2, data = data_reg)
simplex_res_3 <- simplexreg(formula = frl_3, data = data_reg)


## Function to create a coef table for all three regressions
coef_table <- function(tobit_res, beta_res, simplex_res, n_dig = 6){
  ## Coefs
  coef_tobit   <- coeftest(tobit_res)
  coef_beta    <- coeftest(beta_res)
  coef_simplex <- coeftest(simplex_res)
  
  ## Remove phi and log(scale)
  coef_tobit <- coef_tobit[-nrow(coef_tobit),]
  coef_beta  <- coef_beta[-nrow(coef_beta),]
  
  ## format table
  table_tobit <- data.frame(
    Variables = rownames(coef_tobit),
    Tobit = paste(round(coef_tobit[,1], n_dig), " (",
                  round(coef_tobit[,2], n_dig), ") ",
                  ifelse(coef_tobit[,4] < 1e-3, "***",
                         ifelse(coef_tobit[,4] < 1e-2, "**",
                                ifelse(coef_tobit[,4] < 0.05, "*", ""))),
                  sep = "")
  )
  
  table_beta <- data.frame(
    Variables = rownames(coef_beta),
    Beta = paste(round(coef_beta[,1], n_dig), " (",
                 round(coef_beta[,2], n_dig), ") ",
                 ifelse(coef_beta[,4] < 1e-3, "***",
                        ifelse(coef_beta[,4] < 1e-2, "**",
                               ifelse(coef_beta[,4] < 0.05, "*", ""))),
                 sep = "")
  )
  
  table_simplex <- data.frame(
    Variables = rownames(coef_simplex),
    Simplex = paste(round(coef_simplex[,1], n_dig), " (",
                    round(coef_simplex[,2], n_dig), ") ",
                    ifelse(coef_simplex[,4] < 1e-3, "***",
                           ifelse(coef_simplex[,4] < 1e-2, "**",
                                  ifelse(coef_simplex[,4] < 0.05, "*", ""))),
                    sep = "")
  )
  
  ## merge table
  table_final <- merge(table_tobit, table_beta, by = "Variables", 
                       all = TRUE, sort = FALSE)
  table_final <- merge(table_final, table_simplex, by = "Variables", 
                       all = TRUE, sort = FALSE)
  ## return
  return(data.frame(table_final))
  
}

## Table for three regressions
table_res_1 <- coef_table(tobit_res = tobit_res_1,
                          beta_res = beta_res_1,
                          simplex_res = simplex_res_1,
                          n_dig = 2)


table_res_2 <- coef_table(tobit_res = tobit_res_2,
                          beta_res = beta_res_2,
                          simplex_res = simplex_res_2,
                          n_dig = 6)


table_res_3 <- coef_table(tobit_res = tobit_res_3,
                          beta_res = beta_res_3,
                          simplex_res = simplex_res_3,
                          n_dig = 6)

## Format Variables names
table_res_2$Variables <- gsub(pattern = "`", 
                              replacement = "", 
                              x = table_res_2$Variables)
table_res_3$Variables <- gsub(pattern = "`", 
                              replacement = "", 
                              x = table_res_3$Variables)

## Save tables
write_xlsx(x = table_res_1, path = "Robustness Check Regression 1.xlsx")
write_xlsx(x = table_res_2, path = "Robustness Check Regression 2.xlsx")
write_xlsx(x = table_res_3, path = "Robustness Check Regression 3.xlsx")


## Heteroskedasticity check for Robustness Check
plot_robust_reg_check <- function(tobit_res, beta_res, simplex_res){
  require(ggplot2)
  require(reshape2)
  require(ggpubr)
  
  ## Residuals plot
  df_residual <- data.frame(Obs     = 1:length(resid(tobit_res)),
                            Tobit   = resid(tobit_res),
                            Beta    = resid(beta_res),
                            Simplex = resid(simplex_res))
  
  df_residual <- melt(df_residual, id.vars = "Obs")
  
  residual_plot <- ggplot(data = df_residual) + 
    facet_wrap("variable", ncol = 3, scales = "free") +
    geom_point(aes(x = Obs, y = value)) + 
    geom_smooth(aes(x = Obs, y = value), method = "gam") + 
    ylab("Resid") + xlab("Observation") +  
    ggtitle("Resid Plot") + 
    theme_bw()
  
  ## QQ plot
  df_qq <- data.frame(Tobit   = resid(tobit_res),
                      Beta    = resid(beta_res),
                      Simplex = resid(simplex_res))
  df_qq <- melt(df_qq)
  
  plot_qq <- 
    ggplot(mapping = aes(sample = value), data = df_qq) + 
    facet_wrap("variable", ncol = 3, scales = "free") + 
    stat_qq() + 
    stat_qq_line() +
    ylab("Sample") + xlab("Theoretical") + 
    ggtitle("QQ plot") + 
    theme_bw()
  
  ## Check KL divergence of residuals with a normal sample
  KL_residuals <- function(reg_resid){
    ## KL matrix
    kl_mat <- density(reg_resid, n = 1024)
    kl_mat <- data.frame(Quantiles = kl_mat$x,
                         Resid = kl_mat$y,
                         Normal = dnorm(x = kl_mat$x, 
                                        mean = 0, 
                                        sd = sd(reg_resid)))
    
    ## density plot
    kl_gg <- melt(kl_mat, id.vars = "Quantiles")
    Density_Plot <- 
      ggplot(data = kl_gg) +
      geom_line(aes(x = Quantiles, y = value, color = variable)) +
      ylab("Density") + xlab("Quantiles") +
      scale_color_discrete("") + 
      theme_bw()
    
    ## return list
    res <- list()
    res$DensityPlot <- Density_Plot
    res$Divergence <- KLdiv(as.matrix(kl_mat[,-1]))[1,2]
    
    ## return
    return(res)
  }
  
  kl_div_tobit   <- KL_residuals(reg_resid = resid(tobit_res))
  kl_div_beta    <- KL_residuals(reg_resid = resid(beta_res))
  kl_div_simplex <- KL_residuals(reg_resid = resid(simplex_res))
  
  
  ## Density plot of residuals
  density_plot_tobit <- kl_div_tobit$DensityPlot + 
    ggtitle(paste("Tobit Reg. - KL Divergence: ", 
                  round(kl_div_tobit$Divergence, digits = 2),
                  sep = ""))
  density_plot_beta <- kl_div_beta$DensityPlot + 
    ggtitle(paste("Beta Reg. - KL Divergence: ", 
                  round(kl_div_beta$Divergence, digits = 2),
                  sep = ""))
  density_plot_simplex <- kl_div_simplex$DensityPlot + 
    ggtitle(paste("Simplex Reg. - KL Divergence: ", 
                  round(kl_div_simplex$Divergence, digits = 2),
                  sep = "")) 
  
  ## Plot
  heteroske_plot <- ggarrange(
    plotlist = list(residual_plot,
                    plot_qq,
                    ggarrange(plotlist = list(density_plot_tobit,
                                              density_plot_beta,
                                              density_plot_simplex), 
                              ncol = 3, nrow = 1)), 
    nrow = 3
  )
  
  
  ## return
  return(heteroske_plot)
  
}

pdf(file = "Heteroskedasticity Residuals Check Robust Reg Model 1.pdf", 
    width = 10, height = 8)
plot_robust_reg_check(tobit_res = tobit_res_1,
                      beta_res = beta_res_1,
                      simplex_res = simplex_res_1)
dev.off()

pdf(file = "Heteroskedasticity Residuals Check Robust Reg Model 2.pdf", 
    width = 10, height = 8)
plot_robust_reg_check(tobit_res = tobit_res_2,
                      beta_res = beta_res_2,
                      simplex_res = simplex_res_2)
dev.off()

pdf(file = "Heteroskedasticity Residuals Check Robust Reg Model 3.pdf", 
    width = 10, height = 8)
plot_robust_reg_check(tobit_res = tobit_res_3,
                      beta_res = beta_res_3,
                      simplex_res = simplex_res_3)
dev.off()



################################################################################
## Robustness Checks for Region Tobit Regression
################################################################################
## Beta regression
beta_N_NE    <- betareg(formula = frl_region, data = data_reg_N_NE)
beta_S_SE_CW <- betareg(formula = frl_region, data = data_reg_S_SE_CW)


## Simplex regression
simplex_N_NE    <- simplexreg(formula = frl_region, data = data_reg_N_NE)
simplex_S_SE_CW <- simplexreg(formula = frl_region, data = data_reg_S_SE_CW)


## Table for three regressions
table_res_N_NE <- coef_table(tobit_res = tobit_N_NE_res,
                             beta_res = beta_N_NE,
                             simplex_res = simplex_N_NE,
                             n_dig = 6)


table_res_S_SE_CW <- coef_table(tobit_res = tobit_S_SE_CW_res,
                                beta_res = beta_S_SE_CW,
                                simplex_res = simplex_S_SE_CW,
                                n_dig = 6)

## Format Variables names
table_res_N_NE$Variables <- gsub(pattern = "`", 
                                 replacement = "", 
                                 x = table_res_N_NE$Variables)
table_res_S_SE_CW$Variables <- gsub(pattern = "`", 
                                    replacement = "", 
                                    x = table_res_S_SE_CW$Variables)

## Save tables
write_xlsx(x = table_res_N_NE, 
           path = "Robustness Check Region Regression N NE.xlsx")
write_xlsx(x = table_res_S_SE_CW, 
           path = "Robustness Check Region Regression S SE CW.xlsx")


## Heteroskedasticity check for Robustness Check
pdf(file = "Heteroskedasticity Residuals Check Region Robust Reg N NE.pdf", 
    width = 10, height = 8)
plot_robust_reg_check(tobit_res = tobit_N_NE_res,
                      beta_res = beta_N_NE,
                      simplex_res = simplex_N_NE)
dev.off()

pdf(file = "Heteroskedasticity Residuals Check Region Robust Reg S SE CW.pdf", 
    width = 10, height = 8)
plot_robust_reg_check(tobit_res = tobit_S_SE_CW_res,
                      beta_res = beta_S_SE_CW,
                      simplex_res = simplex_S_SE_CW)
dev.off()


################################################################################
## Check natural clustering with PCA
################################################################################
## Principal component analysis
d_pca <- d
colnames(d_pca) <- c("Municipal Spending on Education",
                     "State Spending on Education",    
                     "IDEB Early Years",
                     "IDEB Final Years",                           
                     "IDEB High School",
                     "Literacy Rate",                              
                     "School Attendance",
                     "School Attainment")
d_pca <- principal(r = d_pca, nfactors = ncol(d_pca), scores = TRUE)

## plot the eigenvalues
df_eigen <- data.frame(Components = 1:ncol(d),
                       Eigenvalues = d_pca$values,
                       CumulativeVariance = 100*d_pca$Vaccounted[3,])
colnames(df_eigen) <- c("Components", "Eigenvalues", "Cumulative Variance (%)")
df_eigen <- melt(df_eigen, id.vars = "Components")

pdf(file = "PCA Analysis - Number of factors.pdf", width = 8, height = 4)
ggplot(data = df_eigen) + 
  facet_wrap("variable", ncol = 2, scales = "free") + 
  geom_point(aes(x = Components, y = value)) + 
  geom_line(aes(x = Components, y = value)) + 
  geom_vline(xintercept = 5, lty = "dashed") + 
  ylab("") + 
  theme_bw()
dev.off()

## Components Loadings
n_comp <- 5
pca_loadings <- apply(d_pca$loadings[,1:n_comp], 2, function(x){
  x[abs(x) < 0.5] <- 0
  return(x)
  })
pca_loadings <- data.frame(Variables = rownames(pca_loadings),
                           pca_loadings)
df_loadings_heatmap <- melt(pca_loadings, id.vars = "Variables")

## Heatmap with PCA loadings
pdf(file = "PCA Analysis - Heatmap of loadings.pdf", width = 8, height = 6)
ggplot(df_loadings_heatmap) + 
  geom_tile(aes(x = variable, y = Variables, fill = value)) + 
  scale_fill_gradient("", low = "white", high = "red") + 
  ylab("Variables") + xlab("Component") + 
  ggtitle("Only Variables with absolute loadings greater than 0.5 ") + 
  theme_bw()
dev.off()

## Scores matrix
d_decomp_pca <- data.frame(d_pca$scores[,1:n_comp])
colnames(d_decomp_pca) <- c("Students Performance",
                            "Literacy Rate", 
                            "Spending on Education and Culture (State)",
                            "Spending on Education and Culture (Municipality)",
                            "IDEB High School")
d_decomp_pca <- as.matrix(d_decomp_pca)

################################################################################
## Max number of PCA components to maximize Entropy of TOPSIS scores
################################################################################
## Variables importances signal
imp_pca <- c("+", "+", "-", "-", "+")

## Variable with TOPSIS Score for components
TOPSIS_NG_pca <- NULL

## NG Evaluation with 2,3,4 and 5 components
for(comp in 2:n_comp){
  
  ## data frame with all componenets analysed
  y <- data.frame(d_decomp_pca[ ,1:comp])
  
  product_number  <- nrow(y)  ## Number of DMUs
  criteria_number <- ncol(y)  ## Number of variables
  
  ## Matrices with decisions values
  w    <- matrix(0, nrow = product_number, ncol = criteria_number)
  wres <- matrix(0, nrow = product_number, ncol = criteria_number)
  Ind  <- matrix(0, nrow = product_number, ncol = 1)
  
  ## Matrix with input parameters 
  yn <- data.frame(y)
  
  ## Normalization between 0-1
  yn <- apply(yn, 2, function(x){(x - min(x))/(max(x) - min(x))})
  
  ## Loop in DMU to solve LP
  for (i in 1:product_number) {
    w[i,] <- 0
    lprec <- make.lp(0, criteria_number)
    set.objfn(lprec, c(yn[i,])) ## Objective function
    
    ## First Restriction
    add.constraint(lprec, c(w[i,]+1), "=", 1) ## Sum to one
    
    ## Set of restrictions 2
    for (j in 1:(criteria_number-1)) {
      w[i,]    <-0
      w[i,j]   <- 1
      w[i,j+1] <- -1
      add.constraint(lprec,c(w[i,]), ">=", 0)
    }
    
    ## Set of restrictions 3
    for (j in 1:criteria_number) {
      w[i,]  <-0
      w[i,j] <-1
      add.constraint(lprec,c(w[i,]), ">=", 0)
    }
    
    ## Solve linear programming
    lp.control(lprec,sense='max')
    solve(lprec)
    
    ## add results in matrices
    wres[i,] <- get.variables(lprec)
    Ind[i]   <- get.objective(lprec)
  }
  

  ## Optimized Weights
  w_NG_pca <- wres
  
  ## Mean Weights NG
  w_NG_pca <- colMeans(w_NG_pca)
  
  ## Normalize to sum 1.0
  w_NG_pca <- w_NG_pca/sum(w_NG_pca)
  
  ## Correct names
  names(w_NG_pca) <- colnames(d_decomp_pca[,1:comp])
  
  ## TOPSIS evaluation
  TOPSIS_NG_pca <- cbind(TOPSIS_NG_pca,
                         topsis(decision = d_decomp_pca[,1:comp], 
                                weights = w_NG_pca, 
                                impacts = imp_pca[1:comp])$score)
}

## add original NG
TOPSIS_NG_pca <- cbind(Original = topsis_NG$score,
                       TOPSIS_NG_pca)

## Correct names
colnames(TOPSIS_NG_pca) <- c("Original", 
                             paste(2:n_comp, "Components", sep = " "))

## entropy for all TOPSIS scores
Entropy_PCA_Analysis <- apply(TOPSIS_NG_pca, 2, entropy)

## plot entropy
df_entropy_pca <- melt(Entropy_PCA_Analysis)
df_entropy_pca <- data.frame(Decomposition = rownames(df_entropy_pca),
                             Entropy = df_entropy_pca$value)

pdf(file = "PCA Analysis - Entropy.pdf", width = 6, height = 5)
ggplot(data = df_entropy_pca) + 
  geom_col(aes(x = Decomposition, y = Entropy, fill = Entropy == max(Entropy))) +
  geom_text(aes(x = Decomposition, y = Entropy + 0.2, 
                label = round(Entropy, digits = 3))) +
  theme_bw() + 
  theme(legend.position = "none")
dev.off()


################################################################################
## Descriptive Statistics
################################################################################
descriptive_stats <- do.call(rbind,
                             apply(rank[,c(3:10,16:19)], 2, 
                                   function(x){
                                     data.frame(Min = min(x),
                                                Max = max(x),
                                                Mean = mean(x),
                                                SD = sd(x))
                                   }))

descriptive_stats <- data.frame(Variables = rownames(descriptive_stats),
                                descriptive_stats)


write_xlsx(x = descriptive_stats, path = "Variables Desc Stats.xlsx")

################################################################################
## Appendix table A1 - TOPSIS scores
################################################################################
appendix_table <- data.frame(State = rank[,1], Year = rank[,2],
                             TOPSIS_Entropy = topsis_entropy$score,
                             TOPSIS_Entropy_rank = topsis_entropy$rank,
                             TOPSIS_Ng = topsis_NG$score,
                             TOPSIS_Ng_rank = topsis_NG$rank)

## order with Ng rank
appendix_table <- appendix_table[order(appendix_table$TOPSIS_Ng_rank),]

## Names
colnames(appendix_table) <- c("State", "Year", 
                              "TOPSIS Score (Entropy)", "TOPSIS Rank (Entropy)",
                              "TOPSIS Score (Ng)", "TOPSIS Rank (Ng)")

## Save Appendix Table
write_xlsx(x = appendix_table, path = "Appendix Table A1.xlsx")





