############# Alpha Diversity Wrapper Functions ###########


#function to calculate richness, shannon diversity, and inverse simpson 
#diversity from the phyloseq object and return it all with the sample data as a 
#dataframe

myalphamets<-function(p){ #make a function with one input of a phyloseq object
  
  asvnew<-as.matrix((data.frame(otu_table(p))))#convert otu table to matrix
  
  samdf<-data.frame(sample_data(p))#convert sample data to a dataframe
  
  rich1 <- specnumber(asvnew)#calculates richness
  
  shannon1 <- vegan::diversity(asvnew,"shannon")#calculates the shannon diversity using the function in vegan
  
  invsimp1 <- vegan::diversity(asvnew, "inv")#calculates the inverse simpson diversity using the function in vegan
  
  alp <- cbind(rich1, shannon1, invsimp1,samdf)#put the calculated alpha diversity into the dataframe
  alp <- rownames_to_column(alp, var = "id")
  alp <- data.frame(alp)
  return(alp)
}

### summary statistics for alpha diversity metrics

alpha.summary<- function(ps,sam.var){
  
  alphamets<-myalphamets(ps)
  
  rich.mean<-aggregate(alphamets$rich1~sam.var,FUN=mean)[,2] #mean
  rich.se<-aggregate(alphamets$rich1~sam.var,FUN=function(x) sd(x)/sqrt(length(x)))[,2] #se
  
  shandiv.mean<-aggregate(alphamets$shannon1~sam.var,FUN=mean)[,2] #mean
  shandiv.se<-aggregate(alphamets$shannon1~sam.var,FUN=function(x) sd(x)/sqrt(length(x)))[,2] #se
  
  invsimp.mean<-aggregate(alphamets$invsimp1~sam.var,FUN=mean)[,2] #mean
  invsimp.se<-aggregate(alphamets$invsimp1~sam.var,FUN=function(x) sd(x)/sqrt(length(x)))[,2] #se
  
  all.summary<-data.frame(cbind(levels(as.factor(sam.var)),rich.mean,rich.se,shandiv.mean,shandiv.se,invsimp.mean,invsimp.se))
  return(all.summary)
}

#Samara's function to get summary statistics of the diversity metrics
avgalpha <- function(alpha){ #make a function with 1 input of an alpha diversity dataframe made earlier
  avg<-data.frame(matrix(nrow = 4,ncol = 6)) #make an empty dataframe of 4 rows and 6 columns
  avg[1, ]<-c("alpha diversity metric","Mean","Median","range","standard deviation","standard error") #make a row of labels for the dataframe
  avg[ ,1]<-c("alpha diversity metric","richness","shannon","inv simp") #make a column of labels for the dataframe
  avg[2,2]<-mean(alpha$rich1) #input mean of richness
  avg[2,3]<-median(alpha$rich1) #input median of richness
  range<-data.frame(range(alpha$rich1)) #create a dataframe of range of richness (the range function has 2 outputs)
  range<-as.character(range) #make the values in the range dataframe characters
  avg[2,4]<- as.character(interaction(range,sep="_")) #input the range of richness into the main dataframe merged together
  avg[2,5]<-sd(alpha$rich1) #input the standard deviation of richness
  avg[2,6]<-std.error(alpha$rich1) #input the standard error of richness
  
  avg[3,2]<-mean(alpha$shannon1) #input mean of shannon diversity
  avg[3,3]<-median(alpha$shannon1) #input median of shannon diversity
  range<-data.frame(range(alpha$shannon1)) #create a dataframe of range of shannon diversity
  range<-as.character(range) #make the values in the range dataframe characters
  avg[3,4]<- as.character(interaction(range,sep="_")) #input the range of shannon diversity into the main dataframe merged together
  avg[3,5]<-sd(alpha$shannon1) #input the standard deviation of shannon diversity
  avg[3,6]<-std.error(alpha$shannon1) #input the standard error of shannon diversity
  
  avg[4,2]<-mean(alpha$invsimp1) #input mean of inverse simpson diversity
  avg[4,3]<-median(alpha$invsimp1) #input median of inverse simpson diversity
  range<-data.frame(range(alpha$invsimp1)) #create a dataframe of range of inverse simpson diversity
  range<-as.character(range) #make the values in the range dataframe characters
  avg[4,4]<- as.character(interaction(range,sep="_")) #input the range of inverse simpson diversity into the main dataframe merged together
  avg[4,5]<-sd(alpha$invsimp1) #input the standard deviation of inverse simpson diversity
  avg[4,6]<-std.error(alpha$invsimp1) # input the standard error of inverse simpson diversity
  avg #output the filled dataframe
}
