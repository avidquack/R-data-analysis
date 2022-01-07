setwd("~/OneDrive - Harvard University/Math 156/Datasets/Datasets_export")
setwd("C:/Users/Motol/Desktop/Math 156/datasets")
#Topic 1: Logistic Regression
#Part (a)
#We investigate the relationship between heart disease and age in a random sample of 918 patients who experienced some sort of chest pain
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(resampledata)
heart <- read.csv("heart.csv");head(heart)
age <- heart$Age; disease <- heart$HeartDisease

#disease is a Bernoulli random variable
plot(age, disease)

#Make a contingency table for Age and Heart Disease
tbl<-table(age, disease);head(tbl)

#Turn the vector of ages from the table into a numeric vector
x<-as.numeric(rownames(tbl));x #a vector of the ages

#Now make a table of the probability of getting heart disease as a function of age
prob<-numeric(length(x))
for (i in 1:length(x)) {
  prob[i] = sum((age == x[i])*disease)/sum((age == x[i]))
}
plot(x,prob) 

#Let's find the logistic equation modelling the log-odds of patient getting heart disease against age
#MLE method
library(stats4) 
MLL<- function(a, b)
  -sum(log(exp(a+b*age)/(1+exp(a+b*age)) )*disease+ log(1/(1+exp(a+b*age)))*(1-disease) )
results<-mle(MLL,start = list(a = -0.1, b = -0.02))
results@coef
curve(exp(results@coef[1]+results@coef[2]*x)/ (1+exp(results@coef[1]+results@coef[2]*x)),col = "blue", add=TRUE)

#glm method (yields the same result)
result <-glm(HeartDisease ~ Age, data = heart, family = binomial); result
cfs <- coef(result)
a <- cfs[1]; a
b <- cfs[2]; b
x <- seq(28, 77 , length = 500)
diseaseprob <- function(x,a,b){
  y <- exp(a+b*x)
  y/(1+y)
}
plot(x, diseaseprob(x), ylim = c(0,1), pch = 20, cex = 0.1, xlab = "Age", ylab = "Probability of Heart Disease")

#So for good measure, let's compare the odds of a 30-year old vs 40-year old getting heart disease
p30 <- diseaseprob(30);p30
odds30 <- p30/(1-p30);odds30

p40 <- diseaseprob(40);p40
odds40 <- p40/(1-p40);odds40

ratio <- odds30/odds40;ratio #The odds of a 30 year old getting heart disease is 0.53 times that of a 40 year old ( half as likely!)

#Bootstrap the slope 
N <- 1000
n <- nrow(heart); n
C <- numeric(N)
B <- numeric(N)
A <- numeric(N)
p69 <- numeric(N)
for(i in 1:N){
  index <- sample(n, replace = TRUE)
  heart.boot <- heart[index,]
  result <- glm(HeartDisease ~ Age, data = heart.boot, family = binomial)
  lt <- coef(result)
  A[i] <- lt[1]
  B[i] <- lt[2]
  C[i] <- cor(heart.boot$Age, heart.boot$HeartDisease)
  p69[i] <- diseaseprob(69, lt[1], lt[2])
}
ctrue <- cor(heart.boot$Age, heart.boot$HeartDisease)

hist(C) #histogram of the correlation
abline(v = ctrue)
hist(A) #histogram of the intercept
abline(v = a)
hist(B) #histogram of the slope
abline(v = b, col = "red")

quantile(C, c(.025, .975)) #95% confidence interval for correlations
quantile(A, c(0.025,0.975)) #95% confidence interval for the intercept
quantile(B, c(0.025,0.975)) #95% confidence interval for the slope

p69[i] <- diseaseprob(69, lt[1], lt[2])

#Now let's estimate the probability of a 69 year old with some sort of chest pain having a heart disease
mean(p69) 
q69 <- quantile(p69, probs = c(0.025, 0.925)); q69 #95th percentile interval



#Topic 2: Relative risk (cholesterol and heart disease)

#Exploratory Data Analysis
summary(heart$Cholesterol)
boxplot(heart$Cholesterol)

#Count categories:
nrow(subset(heart, Cholesterol >= 200)) #600 ppl with "high" cholesterol
nrow(subset(heart, Cholesterol >= 0 & Cholesterol < 200)) #318 ppl with "normal" cholesterol

nrow(subset(heart, Cholesterol >= 200 & HeartDisease == 1)) 
nrow(subset(heart, Cholesterol >= 0 & Cholesterol < 200 & HeartDisease == 1)) 

#299 of 600 ppl with high cholesterol got heart disease
#209 of 318 ppl with normal cholesterol got heart disease

h <- 600
n <- 318
p1 <- 299/600;p1
p2 <- 209/318;p2
theta <- p1/p2;theta #0.76, so high cholesterol is less likely than low cholesterol, which seems suspicious


#create data frame
chollevel <- c(rep("high", 600),rep("normal",318))
outcome <- c(rep("hd",299),rep("no hd",301),rep("hd",209),rep("no hd",109))
data <-data.frame(chollevel, outcome);head(data)
table(data$chollevel, data$outcome)

#bootstrap theta (relative risk)
high <- subset(data, select = outcome, subset = (chollevel == "high"), drop = TRUE)
normal <- subset(data, select = outcome, subset = (chollevel == "normal"), drop = TRUE)


#Now do the bootstrapping.
N <- 10^4; ratio <- numeric(N); prop1 <-numeric(N); prop2 <- numeric(N)
for (i in 1:N) {
  sample1 <-sample(high,h, replace = TRUE)
  sample2 <-sample(normal,n, replace = TRUE)
  prop1[i] <- mean(sample1 == "hd")
  prop2[i] <- mean(sample2 == "hd")
  ratio[i] <-prop1[i]/prop2[i]
}
hist(ratio, xlab = "Relative Risk")
abline(v = mean(ratio), col = "red")
abline(v = theta, col = "blue")
bias <- mean(ratio) - theta; bias #0.001317
SE <- sd(ratio); SE #0.044

bias/SE #0.0298



#Topic 3: Chi-Square: Gender and Heart Disease
#We now investigate the independence (chi-square) between heart disease and gender in a random sample of 918 patients

Gender <- heart$Sex
Sick <- heart$HeartDisease
tbl <- table(Gender, Sick);tbl
#Built-in
chisq.test(tbl)

chisq <-function(Obs){
  Expected <- outer(rowSums(Obs), colSums(Obs)/sum(Obs))
  sum((Obs-Expected)^2/Expected)
}
observed <- chisq(tbl); observed

#Create a dataframe
heartdata <- data.frame(Gender, Sick);head(heartdata)

#Permute the gender column
N <- 10^4-1; result <- numeric(N)
for (i in 1:N) {
  Gender.perm<-sample(heartdata$Gender)
  result[i]<-chisq(table(Gender.perm,heartdata$Sick))
}
hist(result, main = "Gender")
abline(v = observed)
pval <- (sum(result >= observed)+1)/(N+1);pval

#We conclude that gender and getting heart disease is not independent (H_a)

#Now to test the relationship between Exercise Angina and Heart disease
Angina = heart$ExerciseAngina
tbl2 <- table(Angina, Sick);tbl2
#Built-in
chisq.test(tbl2)

observed <- chisq(tbl2); observed

#Create a dataframe
heartdata2 <- data.frame(Angina, Sick);head(heartdata2)

#Permute the Exercise Angina column
N <- 10^4-1; result <- numeric(N)
for (i in 1:N) {
  Angina.perm<-sample(Angina, length(Angina))
  result[i]<-chisq(table(Angina.perm,heartdata2$Sick))
}
hist(result, main = "Execise Angina")
abline(v = observed)
pval <- (sum(result >= observed)+1)/(N+1);pval



#Another Chi square test, but now on the relationship between Resting ECG and Heart disease
ECG = heart$RestingECG
tbl3 <- table(ECG, Sick);tbl3
#Built-in
chisq.test(tbl3)

observed <- chisq(tbl3); observed

#Create a dataframe
heartdata3 <- data.frame(ECG, Sick);head(heartdata3)

#Permute the Resting ECG column
N <- 10^4-1; result <- numeric(N)
for (i in 1:N) {
  ECG.perm<-sample(ECG, length(Angina))
  result[i]<-chisq(table(ECG.perm,heartdata2$Sick))
}
hist(result, main = "Resting ECG")
abline(v = observed)
pval <- (sum(result >= observed)+1)/(N+1);pval



#A final Chi-square test on the relatioship between type of chest pain and Heart disease
pain = heart$ChestPainType
tbl4 <- table(pain, Sick);tbl4
#Built-in
chisq.test(tbl4)

observed <- chisq(tbl4); observed

#Create a dataframe
heartdata4 <- data.frame(pain, Sick);head(heartdata4)

#Permute the Resting ECG column
N <- 10^4-1; result <- numeric(N)
for (i in 1:N) {
  pain.perm<-sample(pain, length(Angina))
  result[i]<-chisq(table(pain.perm,heartdata2$Sick))
}
hist(result, main = "Chest pain type")
abline(v = observed)
pval <- (sum(result >= observed)+1)/(N+1);pval

#Topic 4: Multiple linear regression
#We can look at a 


dset <- heart;dset
predictor <- c(4,5)
#Find the regression line by projection
A <- as.matrix(cbind(rep(1,nrow(dset)),dset[,predictor]))
B <- t(A)%*%A
P <- A%*%solve(B)%*%t(A);
y.hat <- P%*%dset$MaxHR ;
lenResid <- sqrt(sum((dset$MaxHR-y.hat)^2));lenResid


plot(dset[,predictor],dset$MaxHR,pch = ".",cex = 3)
#  points(dset[,predictor],y.hat,type = "b")
  
coeff <- solve(B)%*%t(A)%*%dset$MaxHR;coeff

round(coeff[1],3)
round(coeff[1+which.max(predictor==4)],3)
round(coeff[1+which.max(predictor==5)],3)

print("The Equation is", round(coeff[1],3) )

#Topic 5: Bayesian 

nrow(subset(heart, Sex == "M")) #725 men
nrow(subset(heart, Sex == "M" & HeartDisease == 1)) #458

nrow(subset(heart, Sex == "F")) #193 women
nrow(subset(heart, Sex == "F" & HeartDisease == 1)) #50 women

#We assume a flat prior ~ Beta(1,1)

thetam <- 458/725;thetam
thetaw <- 50/193;thetaw

#The likelihood function is 
lm <- thetam^458*(1-thetam)^(725-458);lm
lw <- thetaw^50*(1-thetaw)^(193-50);lw

curve(dbeta(x,1,1))

#The posterior for women is ~ Beta(51, 144)
curve(dbeta(x,51,144))

#So the posterior for men is ~ Beta(459, 277)
curve(dbeta(x,459,277), add = TRUE) 

