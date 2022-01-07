setwd("~/OneDrive - Harvard University/Math 156/Datasets/Datasets_export")
data <- read.csv("covid_datset - covid_datset.csv")

#Topic 1: Contingency Table/Chi-Square independence test

#H0: Party and Vaccination rate is independent
#HA: Party and Vaccination rate is not independent

vcount <- data$tot_pop*0.01*data$fullyv.pop_pct;vcount #calculate total vaccinated
ucount <- data$tot_pop-vcount;ucount #calculate total unvaccinated

data1 <- data.frame(data$Party, vcount, ucount);data1
o11 <- sum(subset(data1, select = "vcount", subset = data.Party == "trump"));o11
o21 <- sum(subset(data1, select = "vcount", subset = data.Party == "biden"));o21
o12 <- sum(subset(data1, select = "ucount", subset = data.Party == "trump"));o12
o22 <- sum(subset(data1, select = "ucount", subset = data.Party == "biden"));o22

table <- matrix(c(o11, o12, o21, o22), ncol=2, byrow=TRUE)
colnames(table) <- c('people vaccinated','people unvaccinated')
rownames(table) <- c('Trump','Biden')
table <- as.table(table)

#Built-in Chi-Square test
chisq.test(table) #p-value 2.2e-16, so we reject H0!

#Topic 2: permutation test by generating subsets by random sampling 

tapply(data$fullyv.pop_pct, data$Party, mean) #Biden states are more vaccinated

index<-which(data$Party=="trump")
observed<-mean(data$fullyv.pop_pct[-index])-mean(data$fullyv.pop_pct[index]); observed
N=10^5-1; result<-numeric(N)  #99999 random subsets
for (i in 1:N) {
  scramble = sample(data$fullyv.pop_pct, length(data$fullyv.pop_pct), replace = FALSE) #random subset
  vaccscram = data.frame(data$Party,scramble)
  result[i]=mean(vaccscram$scramble[-index])-mean(vaccscram$scramble[index])
}

hist(result, breaks = "FD", xlim = range(-15,15))
abline(v = observed, col = "red")
#The observed difference appears extremely unlikely to arise by chance
morevaxxed <-which(result >= observed)
pVal<-(length(morevaxxed)+1)/(length(result)+1) #include the actual sample
pVal 
#EXTREMELY LOW!

#Topic 3: Permutation test of means by creating random subsets

#Let's try 16 states, 8 of which voted for Biden
#Extract the rows and the two columns of interest.
vaxxed <-data[c(1:16), c("Party","fullyv.pop_pct")]; vaxxed 
index <-which(vaxxed$Party=="trump"); index  #use the smaller subset of 8
#Compare mean vaccination rate for the Biden and Trump subsets
observed <-mean(vaxxed$fullyv.pop_pct[-index])-mean(vaxxed$fullyv.pop_pct[index]); observed
#So the states that voted for Biden has larger vaccination rates
#Is it significant?
AllSubsets<-combn(1:16,8)   #there were 8 Trump states
N <-ncol(AllSubsets); N ; diff<-numeric(N) 

#The permutation test
for (i in 1:N) {
  index<-AllSubsets[,i]
  diff[i]<-mean(vaxxed$fullyv.pop_pct[-index])-mean(vaxxed$fullyv.pop_pct[index])
}
hist(diff, breaks = "FD")
abline(v = observed, col = "red")

#The observed gain appears highly unlikely to arise by chance 
better<-which(diff >= observed)   #random subsets that matched or exceeded the observed benefit
pVal<-length(better)/length(diff)
pVal #probability that Biden states have higher vax rates can arise by chance
#EXTREMELY LOW!


#Topic 4: Barplot of region vs. percent fully vaccinated 
table(data$Party, data$Cultural.Region)
barplot(table(data$Party, data$Cultural.Region))


#Topic 5: Histogram of frequency of cases
hist(data$tot_cases, main = "Frequency of total cases", xlab = "Total # of Cases")
hist(data$fullyv.pop_pct, main = "Frequency of Vaccination Rates", xlab = "% Fully Vaccinated")
#Now convert the histogram of frequencies to densities to overlay a normal density function
mu <- mean(data$fullyv.pop_pct);mu
sigma <- sqrt(var(data$fullyv.pop_pct));sigma
hist(data$fullyv.pop_pct, main = "Frequency of Vaccination Rates", xlim = c(0,100), xlab = "% Fully Vaccinated", freq = FALSE)
curve(dnorm(x, mu, sigma), xlim = c(0,100), add= TRUE)

#Sources:
#Populations: https://worldpopulationreview.com/states
#https://data.cdc.gov/browse
#https://www.cnn.com/election/2020/results/president