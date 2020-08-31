install.packages("pROC")
install.packages("dplyr")
install.packages("splitstackshape")
install.packages("moments")
install.packages("sn")

library(rpart)
library(randomForest)
library(pROC)
library(dplyr)
library(splitstackshape) 
library(moments)
library(sn)

set.seed(123)

# Download the data and do some cleaning
df <- read.csv("bank-additional-full.csv", sep=";")
df[df == "unknown"] <- NA
df <- na.omit(df)
df$duration <- NULL # information related to the last call (the dependent variable in the setting)
df$y[df$y == "yes"] <- 1
df$y[df$y == "no"] <- 0
df$y <- as.numeric(df$y)
df$pdays[df$pdays == 999] <- 0 # 999 means that client was not previously contacted
vector.quali <- c('job', 'marital', 'education', 'default', 'housing', 'loan', 'contact',
                  'month', 'day_of_week', 'poutcome', 'y')
df[vector.quali] <- lapply(df[vector.quali], factor)
df[c("age", "campaign", "pdays", "previous")] <- sapply(df[c("age", "campaign", "pdays", "previous")],as.numeric)

# Get fraction (in percentage) of each categorical variable
x <- lapply(df[, c('job', 'marital', 'education', 'default', 'housing', 'loan', 'contact',
                   'month', 'day_of_week', 'poutcome', 'y')], table)

neat.table <- function(x, name){
  xx <- data.frame(x)
  names(xx) <- c("Value", "Count")
  xx$Fraction <- with(xx, Count*100/sum(Count))
  data.frame(Variable = name, xx)
}

do.call(rbind, lapply(seq_along(x), function(i)neat.table(x[i], names(x[i]))))

# information on "time-invariant" continuous variables: age, campaign, pdays, previous
hist(df$previous)
hist(df$age)
hist(df$campaign)
hist(df$pdays)

"""
Fuction for the data simulation.
 
 input:  N (integer): total number of observations to be simulated; should be factor of 77.
         dataset (dataframe): original data set.
         prob.success (numerical): probability threshold for the response variable determening success/failure.
         
 output: dataframe: simulated data set with additional columns of y.prob.true (true probability of success), 
                    y.true (assigned success, i.e 1, if y.prob.true > 0.5), e (error rate), y.prob (sum of y.prob.true and error rate).
 
"""

df.simulation <- function(N, dataset, prob.success=0.5){
  
  # "time-invariant" variables
  # categorical:
  job <- sample(c("housemaid", "services", "admin.", "technician", "blue-collar", "unemployed", "retired", 
                  "entrepreneur", "management", "student", "self-employed"), N, replace=TRUE, 
                prob=c(0.027, 0.09, 0.287, 0.1795, 0.186, 0.024, 0.0399, 0.0357, 0.0758, 0.02, 0.0358))
  marital <- sample(c("divorced", "married", "single"), N, replace=TRUE, prob=c(0.1165, 0.5737, 0.3097))
  education <- sample(c("basic.4y", "basic.6y", "basic.9y", "high.school", "illiterate", 
                        "professional.course", "university.degree"), N, replace=TRUE, 
                      prob=c(0.078, 0.046, 0.14, 0.25, 0.0004, 0.14, 0.342))
  default <- sample(c("no", "yes"), N, replace=TRUE, prob=c(99.99, 0.01))
  housing <- sample(c("no", "yes"), N, replace=TRUE, prob=c(45.8, 54.2))
  loan <- sample(c("no", "yes"), N, replace=TRUE, prob=c(84.36, 15.64))
  contact <- sample(c("cellular", "telephone"), N, replace=TRUE, prob=c(67.05, 32.95))
  day_of_week <- sample(c("fri", "mon", "thu", "tue", "wed"), N, replace=TRUE, prob=c(0.188, 0.206, 0.209, 
                                                                                      0.195, 0.2))
  poutcome <- sample(c("failure", "nonexistent", "success"), N, replace=TRUE, prob=c(0.114, 0.847, 0.39))
  #numerical:
  age.params <- cp2dp(c(mean(dataset$age), sd(dataset$age), skewness(dataset$age)), "SN")
  age <- c(rsn(N, dp=age.params))
  campaign.params <- cp2dp(c(mean(dataset$age), sd(dataset$age), skewness(dataset$age), 
                             kurtosis(dataset$campaign)), "ST")
  campaign <- c(rst(N, dp=campaign.params))
  pdays.params <- cp2dp(c(mean(dataset$pdays), sd(dataset$pdays), 3.3, 
                          kurtosis(dataset$pdays)), "ST")
  pdays <- c(rst(N, dp=pdays.params))
  previous.params <- cp2dp(c(mean(dataset$previous), sd(dataset$previous), 2.5, 
                             kurtosis(dataset$previous)), "ST")
  previous <- c(rst(N, dp=previous.params))
  
  #create dataframe out of "time-invariant" variables
  df.nottime.var <- cbind(job, marital, education, default, housing, loan, contact, day_of_week, 
                          poutcome, age, campaign, pdays, previous)
  colnames(df.nottime.var) <- c("job", "marital", "education", "default", "housing", "loan",
                                "contact", "day_of_week", "poutcome", "age", "campaign", 
                                "pdays", "previous")
  df.nottime.var <- as.data.frame(df.nottime.var)
  
  # time-dependent variables
  # Create a data frame which contains time-dependent variables. 
  # The values are chosen randomly  for each quarter from the original data set. 
  # Each quarter contains the same number of people interviewed (factor of 7).
  df.time.unique <- unique(dataset[,c('month', 'emp.var.rate','nr.employed','cons.conf.idx','cons.price.idx','euribor3m')])
  df.time.strat <- stratified(stratified(df.time.unique, "cons.conf.idx", 7), "nr.employed", 7)
  n.per.time <- N/nrow(df.time.strat)
  df.time.var <- df.time.strat[rep(seq_len(nrow(df.time.strat)), each = n.per.time), ]
  
  # dataframe of all explanatory variables
  df.full <- cbind(df.nottime.var, df.time.var)
  df.full$age <- as.integer(df.full$age)
  df.full$campaign <- as.integer(df.full$campaign)
  df.full$pdays <- as.integer(df.full$pdays)
  df.full$previous <- as.integer(df.full$previous)
  
  # use the whole empirical dataset to create tree and produce y.true for the simulated data set.
  tree.fit <- rpart(y~., data=dataset, method="class")
  y.predict.matrix <- predict(tree.fit, newdata=df.full, type="prob")
  df.full$y.prob.true <- as.vector(y.predict.matrix[,2])
  df.full$y.true <- as.factor(ifelse(df.full$y.prob.true > 0.5, 1, 0))
  df.full$e <- rnorm(n=N, mean=0, sd=0.2)
  df.full$y.prob <- df.full$y.prob.true + df.full$e
  df.full$y <- as.factor(ifelse(df.full$y.prob > prob.success, 1, 0))
  
  return(df.full)
  
}  
  
  
"""
Function for calculating  test Area Under the Curve (AUC) scores 
of Decision Tree and Random Forest applied to the simulated data set
using a rolling window evaluation. 

input: dataset (dataframe): simulated data set
       K (integer): number of runs, i.e. training and test sets updates
       test.size (numerical): fraction of the data used as a test set; 0 < test.size < 1

output: dataframe: results for the average test Area Under the Curve (AUC)

"""
  sim.rolling.window.evaluation <- function(dataset, K, test.size){
    
    # calculate sizes of training and test sets
    r <- (1-test.size)/test.size
    N <- nrow(dataset) 
    n.train <- N*r/(K+r)
    n.test <- (N-n.train)/ K  
    
    # additional data needed for simulation evaluation
    y.true <- dataset$y.true
    y.true.train <- head(x=y.true, n=n.train)
    y.true.test.all <- tail(x=y.true, n=(N-n.train))
    y.true.test <- head(x=y.true.test.all,n=n.test)
    # drop unnecessary columns
    dataset$y.prob.true <- NULL
    dataset$e <- NULL
    dataset$y.prob.fin <- NULL
    dataset$y.true <- NULL
    
    # create: initial training set (train),
    #         set from which test sets are drawn (test.all), 
    #         first test set  (test)
    
    train <- head(x=dataset, n=n.train)
    y.train <- train$y
    predictors.train <- train[, names(train)!= "y"]
    test.all <- tail(x=dataset, n=(N-n.train))
    test <- head(x=test.all,n=n.test)
    y.test <- test$y
    predictors.test <- test[, names(test)!= "y"]
    
    
    # loop for 1) fitting each model using a train set, calculating respective test AUC 
    #          2) updating train and test data sets
    #          3) calculating test AUC
    
    AUC.tree.appended <- c()
    AUC.forest.appended <- c()
    
    for (i in 2:K){
      
      forest.fit <- randomForest(x=predictors.train, y=y.train)
      y.matrix.forest <- predict(object=forest.fit, newdata=predictors.test, type="prob")
      y.data.forest <- cbind.data.frame(as.vector(y.matrix.forest[,2]), y.true.test)
      AUC.forest.result <- auc(data=y.data.forest, response=y.data[,2], predictor=y.data[,1])
      AUC.forest.appended <- append(AUC.forest.appended, AUC.forest.result)
      
      tree.fit <- rpart(y~., data=train, method="class")
      y.matrix.tree <- predict(tree.fit, newdata=test, type="prob")
      y.data.tree <- cbind.data.frame(as.vector(y.matrix.tree[,2]), y.true.test)
      AUC.tree.result <- auc(data=y.data.tree, response=y.data[,2], predictor=y.data[,1])
      AUC.tree.appended <- append(AUC.tree.appended, AUC.tree.result)
      
      if (nrow(test.all) >= n.test){
        
        train <- tail(train, -n.test)
        train <- rbind(train, test)
        test.all <- tail(test.all, -n.test)
        test <- head(x=test.all,n=n.test)
        y.train <- train$y
        predictors.train <- train[, names(train)!= "y"]
        y.test <- test$y
        predictors.test <- test[, names(test)!= "y"]
        
        y.true.test.all <- tail(y.true.test.all, -n.test)
        y.true.test <- head(x=y.true.test.all,n=n.test)
        
      }
      
    }
    
    ave.AUC.forest <- mean(AUC.forest.appended)
    ave.AUC.tree <- mean(AUC.tree.appended)
    
    df.AUC <- data.frame(matrix(ncol=2, nrow=2))
    colnames(df.AUC) <- c("Method","ave.AUC")
    df.AUC$Method[1] <- "Decision Tree"
    df.AUC$Method[2] <- "Random Forest"
    df.AUC$ave.AUC[1] <- ave.AUC.tree
    df.AUC$ave.AUC[2] <- ave.AUC.forest
    
    return(df.AUC)
    
  }
  
  # Simulate data set
  sim.data <- df.simulation(N=7700, dataset=df)
  
  # Get the evaluation results
  sim.results <- sim.rolling.window.evaluation(df.trial, 10, 0.10)
  print(sim.results)
  
  # remove unnecessary columns for the further analysis
  sim.data$y.prob.true <- NULL
  sim.data$e <- NULL
  sim.data$y.prob.fin <- NULL
  sim.data$y.true <- NULL
  
  # fit Classification tree and get summary
  tree.fit <- rpart(y~., data=sim.data, method="class")
  summary(tree.fit)
  plot(tree.fit)
  text(tree.fit)
  
  #check for multicollinearity of a subset of features
  df.col <- data.frame(sim.data$emp.var.rate, sim.data$euribor3m, sim.data$nr.employed)
  df.col[] <- lapply(df.col, as.numeric)  
  chisq <- chisq.test(table(unlist(df.col)))
  print(chisq)
  
  df.col.add <- data.frame(sim.data$emp.var.rate, sim.data$euribor3m, sim.data$nr.employed, sim.data$cons.conf.idx,  sim.data$cons.price.idx)
  df.col.add[] <- lapply(df.col.add, as.numeric)  
  chisq.add <- chisq.test(table(unlist(df.col.add)))
  print(chisq.add)

  
"""
Function for calculating  test Area Under the Curve (AUC) scores 
of Decision Tree and Random Forest applied to the original data set
using a rolling window evaluation. 

input: dataset (dataframe): data set
       K (integer): number of runs, i.e. training and test sets updates
       test.size (numerical): fraction of the data used as a test set; 0 < test.size < 1

output: dataframe: results for the average test Area Under the Curve (AUC)

"""
  rolling.window.evaluation <- function(dataset, K, test.size){
    
    # calculate size of a training set (n.train)
    r <- (1-test.size)/test.size
    N <- nrow(dataset) 
    n.train <- N*r/(K+r)
    
    
    # create: initial training set (train),
    #         set from which test sets are drawn (test.all), 
    #         first test set  (test)
    
    train <- head(x=dataset, n=n.train)
    y.train <- train$y
    predictors.train <- train[, names(train)!= "y"]
    test.all <- tail(x=dataset, n=(N-n.train))
    n.test <- nrow(test.all)/ K  # size of each test set
    test <- head(x=test.all,n=n.test)
    y.test <- test$y
    predictors.test <- test[, names(test)!= "y"]
    
    
    # loop for 1) fitting each model using a train set, calculating respective test AUC 
    #          2) updating train and test data sets
    #          3) calculating test AUC
    
    AUC.tree.appended <- c()
    AUC.forest.appended <- c()
    
    for (i in 2:K){
      
      forest.fit <- randomForest(x=predictors.train, y=y.train)
      y.matrix.forest <- predict(object=forest.fit, newdata=predictors.test, type="prob")
      y.data.forest <- cbind.data.frame(as.vector(y.matrix.forest[,2]), y.test)
      AUC.forest.result <- auc(data=y.data.forest, response=y.data[,2], predictor=y.data[,1])
      AUC.forest.appended <- append(AUC.forest.appended, AUC.forest.result)
      
      tree.fit <- rpart(y~., data=train, method="class")
      y.matrix.tree <- predict(tree.fit, newdata=test, type="prob")
      y.data.tree <- cbind.data.frame(as.vector(y.matrix.tree[,2]), y.test)
      AUC.tree.result <- auc(data=y.data.tree, response=y.data[,2], predictor=y.data[,1])
      AUC.tree.appended <- append(AUC.tree.appended, AUC.tree.result)
      
      if (nrow(test.all) >= n.test){
        
        train <- tail(train, -n.test)
        train <- rbind(train, test)
        test.all <- tail(test.all, -n.test)
        test <- head(x=test.all,n=n.test)
        y.train <- train$y
        predictors.train <- train[, names(train)!= "y"]
        y.test <- test$y
        predictors.test <- test[, names(test)!= "y"]
        
      }
      
    }
    
    ave.AUC.forest <- mean(AUC.forest.appended)
    ave.AUC.tree <- mean(AUC.tree.appended)
    
    df.AUC <- data.frame(matrix(ncol=2, nrow=2))
    colnames(df.AUC) <- c("Method","ave.AUC")
    df.AUC$Method[1] <- "Decision Tree"
    df.AUC$Method[2] <- "Random Forest"
    df.AUC$ave.AUC[1] <- ave.AUC.tree
    df.AUC$ave.AUC[2] <- ave.AUC.forest
    
    return(df.AUC)
    
  }  
  
  # get rolling window evaluation results
  emp.results <- rolling.window.evaluation(dataset=df,K=10, test.size = 0.10)
  print(emp.results)
  
  # fit Classification tree and get summary
  emp.tree.fit <- rpart(y~., data=df, method="class")
  summary(emp.tree.fit)
  plot(emp.fit.tree)
  text(emp.fit.tree)
  
  # fit Random Forest
  df.response <- df$y
  df$y <- NULL
  forest.fit <- randomForest(x=df, y=df.response)
  importance(forest.fit)
  varImpPlot(forest.fit)