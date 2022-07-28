# These R functions can be used to estimate prevalence from pools
#liklihood function

calc.log.lik<-function(p,hitFlies,missFlies){
        p.misses<-log((1-p)^missFlies)
        p.hits<-log(1-((1-p)^hitFlies))
        return(sum(c(p.hits,p.misses)))
}


########### Plotting and calculating function #########
##takes two numeric vectors, one of 'pot sizes' one of '1|0' for infection status. 
#Plot is optional (set to FALSE if you don't want a pretty picture of the log likelihood profile

calcVbounds<-function(nFlies,hits,plot=FALSE){
        
        hits<-as.numeric(hits[!is.na(hits)])
        nFlies<-as.numeric(nFlies[!is.na(hits)])
        
        #seperate out the pots into hits and misses
        hitFlies<-nFlies[which(as.logical(hits))]
        missFlies<-nFlies[which(!hits)]
        
        #the surface
        results<-unlist(lapply(seq(0,1,0.0001),calc.log.lik,hitFlies,missFlies))
        lower<-min(seq(0,1,0.0001)[which((results-max(results)+2)>0)])
        upper<-max(seq(0,1,0.0001)[which((results-max(results)+2)>0)])
        point<-seq(0,1,0.0001)[which.max(results)]
        
        if(plot==TRUE){
                #plot the surface
                plot(seq(0,1,0.0001),results,type="l",ylim=c((max(results)-10),max(results)),xlab="Prevelance",ylab="log Likelihood")
                abline(h=(max(results)-2))
                abline(v=lower,col="red",lty=3)
                abline(v=upper,col="red",lty=3)
                abline(v=point,col="red",lwd=4)
        }
        return(c(lower,point,upper))
        
}

#Then, after setting up these functions, you can estimate prevalence with

pool_sizes<-c(30, 7, 30, 7, 30, 7, 30, 7, 30, 7, 30, 7, 30, 7, 15)
positive_pools<-c(1,0,1,0,1,0,2,0,1,0,0,0,0,0,1)
calcVbounds(pool_sizes, positive_pools,plot=TRUE)


pool_sizes<-c(30, 7)
positive_pools<-c(2,0)
calcVbounds(pool_sizes, positive_pools,plot=TRUE)


#As an output it gives you {lower bound, ML estimate, upper bound}

# Load the library 
install.packages("epiR")
library(epiR)
library(ggplot2)
library(scales)
library(zoo)

ncas <- 4; npop <- 200
tmp <- as.matrix(cbind(ncas, npop))
epi.conf(tmp, ctype = "prevalence", method = "exact", N = 1000, design = 1, 
         conf.level = 0.95) * 100


# Determining the prevalence 
## Install and load the require packages for the analysis

install.packages("rprev")
library(rprev)
library(survival)

# Importing data

library(readr)
droso <- read_csv("data.csv")
View(droso)

# Calculate the prevalence

summary(prevsim)
prevalence_total <- prevalence(index='2013-01-30', 
                               num_years_to_estimate=c(3, 5, 10, 20), 
                               data=prevsim, 
                               inc_formula = entrydate ~ sex,
                               surv_formula = Surv(time, status) ~ age + sex, 
                               dist='weibull', 
                               population_size = 1e6,
                               death_column = 'eventdate')

