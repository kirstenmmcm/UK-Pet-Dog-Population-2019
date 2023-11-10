# IMPORT DATA ####
## Raw UK pet dog dataset: de-duplicated and source anonymized #### 
FINAL_ALL_DEDUP <- read.delim("UKPDP_DEDUP_ANON_RAW_DATA.txt")

cols1 <- c("PROVIDER","BREED","CROSSBREED","SEX","POSTCODE_AREA","STATUS","MICROCHIP")
#cols2 <- c("DATE_OF_BIRTH","TERMINATION_DATE")

library(dplyr)
FINAL_ALL_DEDUP <- FINAL_ALL_DEDUP %>% 
                    mutate_at(cols1, funs(factor(.))) 

### Subset to alive only and under 18.3 years (95% dead in survival analysis - McMillan et al., submitted) ####
FINAL_ALL_DEDUP_sub<-subset(FINAL_ALL_DEDUP, AGE_YRS <= 18.3 & STATUS == "Alive")

## Previous estimate of dogs per postcode (from Aegerter et al. 2017) ####
#https://data.gov.uk/dataset/ec8fc820-2e36-49d0-a09c-e2901e10b2e4/dog-population-per-postcode-district
dog_estimate <- read.csv("Estimated_dog_pop_per_postcode.csv", check.names=FALSE, na.strings=c("","NA"))
N_sum<-sum(as.integer(dog_estimate$TOTAL_dogs)) # 11,884,824

## Human population per postcode ####
human_population <- read.csv("UK human population by postcode.csv", check.names=FALSE, na.strings=c("","NA"))
human_pop_scaled<-as.numeric(scale(human_population[,2]))
N_sum_human<-sum(as.integer(human_population$human_population_2011)) # 63,305,852


# CREATE MATRIX ####
unique(FINAL_ALL_DEDUP$POSTCODE_AREA) #124 postcodes

dog_per_postcode_per_survey<-table(FINAL_ALL_DEDUP$POSTCODE_AREA, FINAL_ALL_DEDUP$PROVIDER) # 11 providers
n_ij_dedup<-as.matrix.data.frame(dog_per_postcode_per_survey)
  write.csv(n_ij_dedup,'n_ij_dedup.csv')

## All Dogs Matrix ####
    all_ij<-read.csv('n_ij_dedup.csv',header=T)
    all_ij<-all_ij[,-1]
     
## Total Count of Dogs for Poisson-Gamma Mixture Models ####  
    Cval=sum(all_ij)
    
## RNG Seed ####
    set.seed(123)
    
## Function for Quick Plotting of Distribution of N ####    
    n_plot<-function(jagsobj){
      
      library(dplyr)
      library(tidyr)
      library(ggplot2)
      library(tidybayes)
      
      jags_mcmc<-as.data.frame(jagsobj$BUGSoutput$sims.matrix)
      
      ## Long Format #### 
      jags_mcmc %>% pivot_longer(colnames(jags_mcmc),names_to="parameter",values_to="value") -> jags_mcmc_plotdata
      
      #Plot ####
      jags_mcmc_plotdata %>% filter(parameter=="N") %>% ggplot(aes(y=parameter,x=value)) + stat_halfeyeh(shape=21,point_fill="white",point_size=5,colour="black") + theme_bw()
      
    }
    
# PACKAGES FOR PLOTTING AND DATA MANIPULATION ####
    library(ggplot2)
    library(cowplot)
    library(dplyr)
    library(tidybayes)
    library(tidyr)
    library(tibble)
    options(scipen = 5) #suppress abbreviated numbers
    library(runjags)
    
    library(scales) #for adding commas to long numbers on axis labels 
    options(mc.cores=8)
    



     
     
# MODEL: N-MIXTURE MODEL WITH UNCORRELATED RANDOM INTERCEPTS AND SLOPES ####
    ## Overdispersed Poisson with Observation Level Random Effect ####
     
     model_nmix_random = '
                    model {
                      for(i in 1:I) { # Loop through postcode
                        for(j in 1:J) { # Loop through survey
                          n_ij[i,j] ~ dbinom(p_ij[i,j], N_i[i])
                          logit(p_ij[i,j]) = alpha[j] + beta[j] * human_pop[i]
                        }
                        N_i[i] ~ dpois(lambda[i]) # lambda is the prior estimate of number of dogs in this postcode
                        log(lambda[i]) <-  alpha.mu[i]
                          }
                      
                    #PostCode Random Intercepts  
                      for(i in 1:I){
                          alpha.mu[i]~dnorm(0,tau.olre)
                        }
                        
                        sd.olre ~ dunif(0,2)
                        tau.olre <- 1/(sd.olre*sd.olre)
                        
                    
                        #beta.mu ~ dnorm(0,0.001)
                        
                    #Survey Random Intercepts and Slopes   
                      for (j in 1:J) {
                        alpha[j] ~ dnorm(0, tau_p_intercept)
                        beta[j] ~ dnorm(0, tau_p_slope)
                      }
                      
                        sd_p_intercept ~ dunif(0,2)
                        sd_p_slope ~ dunif(0,2)
                        tau_p_intercept <- 1/(sd_p_intercept*sd_p_intercept)
                        tau_p_slope <- 1/(sd_p_slope*sd_p_slope)

                      N = sum(N_i) # Sum up all the populations in each postcode to give total estimate
                    }
                    
                           ' ##end model 
     
     
     ## Inits ~~~~  
     Nst <- apply(all_ij, 1, sum) + 1 #maximum count per row of n_ij
     inits_od <- function(){list(N_i = Nst, alpha.mu=rnorm(124, 0, 1), 
                                 alpha=rnorm(ncol(all_ij), 0, 1), beta=rnorm(ncol(all_ij), 0, 1))}
     
     ## Run Model ####
     jags_nmix_random = jags(data = list(I = nrow(all_ij), # The number of post-codes
                                    J = ncol(all_ij), # The number of surveys
                                    human_pop=human_pop_scaled,
                                    n_ij = all_ij), # Counts of dogs in each postcode (row) and each survey (column)
                        inits=inits_od, #starting values
                        parameters.to.save = c('N', 'alpha','beta', 'alpha.mu','beta.mu','lambda',"N_i","sd.olre","sd_p_intercept","sd_p_slope"),
                        model.file = textConnection(model_nmix_random),n.iter=100000,n.burnin=10000,n.thin=50)
     
     
     nmix_random_out<-  jags_nmix_random$BUGSoutput$summary
     write.csv(nmix_random_out, file = "~/nmix_random_out_summary.csv")
     
      head(nmix_random_out)
     
     
     
     
    ## Strip Out MCMC Chains #### 
      nmix_random_mcmc<-as.data.frame(jags_nmix_random$BUGSoutput$sims.matrix)
      head(nmix_random_mcmc)
      
    ## Rename N_i with Postcodes #### 
      colnames(nmix_random_mcmc)[grepl("N_i",colnames(nmix_random_mcmc))]<-rownames(dog_per_postcode_per_survey)
      colnames(nmix_random_mcmc)[1]<-"N_total"
      
    ## Long Format #### 
      nmix_random_mcmc %>% pivot_longer(colnames(nmix_random_mcmc),names_to="parameter",values_to="value") -> jags_nmix_random_plotdata
      
    ## Plot #### 
      jags_nmix_random_plotdata %>% filter(parameter=="N_total") %>% 
        ggplot(aes(y=parameter,x=value)) + 
        stat_dotsinterval() + #stat_halfeye(shape=21,point_fill="white",point_size=5,colour="black") 
        theme_bw()+
        xlab("Predicted Population Size, Per Postcode")+
        ylab("")+
        theme(axis.text.x = element_text(size=8, angle = 45, hjust = 1, colour="black"),
              axis.text.y = element_text(size=8, hjust = 1, colour="black"),
              axis.title.x = element_text(size=8, colour="black"),
              axis.title.y = element_text(size=8, colour="black"))
      
    ### Plot of Postcode Specific #### 
      jags_postcode_nmix_plotdata <- jags_nmix_random_plotdata %>% filter(parameter %in% rownames(dog_per_postcode_per_survey))
      pc1<-ggplot(jags_postcode_nmix_plotdata,aes(y=parameter,x=value)) + stat_halfeye(shape=21,point_fill="white",point_size=5,colour="black") + theme_bw()
      pc1
      
      
    ## Dataframe Assembly ####
      # Model Posterior Means 
      # Raw Input Data 
      # Previous Per-Postcode Estimates
      
    
      
    ### Summarise Model Means ####
      pc_aggregate_data<-jags_postcode_nmix_plotdata %>% group_by(parameter) %>% summarize(posterior_mean= mean(value))
      
    ### Calculate Highest Posterior Density Intervals from MCMC chains (defaults to 95%) ####
      nmix_random_95<-as.data.frame(HPDinterval(as.mcmc(jags_nmix_random$BUGSoutput$sims.matrix)))
        head(nmix_random_95)
        
    ### Add to Dataframe ####
        pc_aggregate_data$lower95<- nmix_random_95[grep("N_",rownames(nmix_random_95)),1]
        pc_aggregate_data$upper95<- nmix_random_95[grep("N_",rownames(nmix_random_95)),2]
        
     ### Add Input Data ####  
      pc_aggregate_data$raw_data<- apply(all_ij,1,sum)  
      head(pc_aggregate_data)
      
    ### Add Previous Estimates ####    
      pc_aggregate_data$previous_estimate_Aegerter<- dog_estimate[,2]  
      
    ### Add Human Pop Size ####      
      pc_aggregate_data $human_pop2011<- human_population[,2]
      head(pc_aggregate_data)
      
      
  # PLOTS ####
      
    ## Our model versus Aegerter - look for Outliers ####
      nmix_aegerter1<- ggplot(pc_aggregate_data,aes(x=previous_estimate_Aegerter,y=posterior_mean))  + theme_bw() + geom_abline(intercept=0,slope=1)
        nmix_aegerter2<- nmix_aegerter1 + geom_errorbar(aes(x=previous_estimate_Aegerter,ymin=lower95,ymax=upper95),width=0.15)
        nmix_aegerter3<- nmix_aegerter2 + geom_point(size=5,shape=21,fill="white") + scale_y_continuous(labels=scales::comma) + scale_x_continuous(labels=scales::comma)
        nmix_aegerter4<- nmix_aegerter3 + labs(x="Dogs Per Postcode (Aegerter)",y="Dogs Per Postcode (N-Mixture Model)") + theme(axis.text = element_text(size=20),axis.title=element_text(size=20))
        nmix_aegerter4
          ggsave2('Dogs Per Postcode Versus Aegerter.pdf',nmix_aegerter4)
          
    ## Our model versus Human Pop ####
      nmix_human1<- ggplot(pc_aggregate_data,aes(x=human_pop2011,y=posterior_mean)) + theme_bw() + geom_smooth(method="lm",formula= y ~ x +I(x^2))
        nmix_human2<- nmix_human1 + geom_errorbar(aes(x=human_pop2011,ymin=lower95,ymax=upper95),width=0.2)
        nmix_human3<- nmix_human2 + geom_point(size=5,shape=21,fill="white")  + scale_y_continuous(labels=scales::comma) + scale_x_continuous(labels=scales::comma)
        nmix_human4<- nmix_human3 +  labs(x="Per-Postcode Human Population Size (2011 Census)",y="Dogs Per Postcode") + theme(axis.text = element_text(size=20),axis.title=element_text(size=16),axis.text.x = element_text(
          angle=45,hjust=1))
        nmix_human4
            ggsave2('Dogs Per Postcode Versus Human Pop Size.pdf',nmix_human4)
     
      
  ## Dogs Per Capita Per Postcode ####
        head(nmix_random_mcmc)
      
      #Pull Out Per-Postcode Raw Chains // Transpose to allow vectorisation with human pop data
        perpost_raw<-t(as.matrix(nmix_random_mcmc[,which(colnames(nmix_random_mcmc) %in% pc_aggregate_data$parameter)]))
          dim(perpost_raw)
      
      #Divide Chains By Human Pop Size 
         dogs_per_capita_raw<- perpost_raw / pc_aggregate_data$human_pop2011
        
      #Calcualate Means and 95% Credible Intervals
         dogs_per_capita95<-apply(dogs_per_capita_raw,1,function(x)quantile(x,c(0.025,0.975)))
         dogs_per_capita_stats<- data.frame(posterior_mean = apply(dogs_per_capita_raw,1,mean),lower95=dogs_per_capita95[1,],upper95=dogs_per_capita95[2,])
          head(dogs_per_capita_stats)
      
      #Add In Human Pop Size 
          dogs_per_capita_stats$human_pop=pc_aggregate_data$human_pop2011
      
      #Plot
        dogs_capita_plot1<- ggplot(dogs_per_capita_stats,aes(x=human_pop,y=posterior_mean))  
          dogs_capita_plot2<- dogs_capita_plot1 +  geom_errorbar(aes(x=human_pop,ymin=lower95,ymax=upper95),width=0.2)
          dogs_capita_plot3<- dogs_capita_plot2 + geom_point(size=5,shape=21,fill="white") + scale_x_continuous(labels=scales::comma)
          dogs_capita_plot4<- dogs_capita_plot3 + labs(x="Per-Postcode Human Population Size (2011 Census)",y="Dogs Per Capita") + theme(axis.text = element_text(size=20),axis.title=element_text(size=16),axis.text.x = element_text(
            angle=45,hjust=1))
          dogs_capita_plot4
            ggsave2('Dogs Per Capita Versus Human Pop Size.pdf',dogs_capita_plot4)