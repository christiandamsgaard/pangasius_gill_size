setwd("/Users/au231308/Documents/Projects/Pangasius AB balance ILCM/")                                    # Set working directory


# Load packages
Packages <- c("RColorBrewer","ggplot2", "gridExtra","readxl","nlme","lme4","multcomp","multcompView", "cowplot","ggsignif")
lapply(X = Packages, FUN = library, character.only = TRUE) # load all required packages

# Colours
cols                     <- brewer.pal(9,"Set1")


## Custom functions
stderr                   <- function(x) sd(x)/sqrt(length(x))                                      # Standard error of mean
lower                    <- function(x) mean(x)-stderr(x)                                          # Mean - standard error of mean
upper                    <- function(x) mean(x)+stderr(x)                                          # Mean + standard error of mean


# Constants
aco2 <-0.0374 # CO2 solubility at 28 degrees in mmol mmHg-1 (Dejours 1981)


# Load data
data<-read_excel(path = "./data_panga_gill.xlsx", sheet = "Ion regulation") # Import data
data<-as.data.frame(data)

# Calculated oxygen and hemoglobin values
data$pk   <- 4.83+0.17*data$ph 
data$pco2<-data$hco3/(aco2*10^(data$ph-data$pk))

# Calculate compensation
data.comp<-expand.grid(unique(data$fishID),c(4,8,18,24),NA,NA)
colnames(data.comp)<-c("fishID","time","treatment","compensation")
data.comp$fishID<-as.character(data.comp$fishID)

for (k in 1:dim(data.comp)[1]){
  # assign treatment
  data.comp$treatment[k]<-ifelse(test = which(unique(data.comp$fishID)==data.comp$fishID[k])<7,yes = "H", no = "L")
  
  # data frame for the fish individual of interest
  temp.df<-data[data$fishID==data.comp$fishID[k],]
  t0<-temp.df$ph[which(temp.df$time==0)]
  t4<-temp.df$ph[which(temp.df$time==4)]
  t8<-temp.df$ph[which(temp.df$time==8)]
  t18<-temp.df$ph[which(temp.df$time==18)]
  t24<-temp.df$ph[which(temp.df$time==24)]
  range<-t0-t4;range  
  # compensation
  if(data.comp$time[k]==4){
    data.comp$compensation[k]<-0
  }
  
  if(data.comp$time[k]==8){
    data.comp$compensation[k]<-(t8-t4)/range
  }
  
  if(data.comp$time[k]==18){
    data.comp$compensation[k]<-(t18-t4)/range
  }
  
  if(data.comp$time[k]==24){
    data.comp$compensation[k]<-(t24-t4)/range
  }
}
data.comp$compensation<-data.comp$compensation*100

# How much is pH depressed after 4 h?
ph.depression<-
  data.frame(
    fish = unique(data$fishID),
    dph = rep(NA,12),
    treat = c(rep("H",6),rep("L",6))
  )

for (i in 1:12){
  temp<-data[data$fishID==unique(data$fishID)[i],]
  ph.depression$dph[i] <- temp$ph[temp$time==0]-temp$ph[temp$time==4]  
}

aggregate(dph~treat,data = ph.depression,mean)
aggregate(dph~treat,data = ph.depression,stderr)
t.test(ph.depression[-9,]$dph~ph.depression[-9,]$treat)


# How much does HCO3 increase from 4 to 24h?
hco3.accu<-
  data.frame(
    fish = unique(data$fishID),
    dhco3 = rep(NA,12),
    treat = c(rep("H",6),rep("L",6))
  )

for (i in 1:12){
  temp<-data[data$fishID==unique(data$fishID)[i],]
  hco3.accu$dhco3[i] <- temp$hco3[temp$time==24]-temp$hco3[temp$time==4]  
}

aggregate(dhco3~treat,data = hco3.accu,mean)
aggregate(dhco3~treat,data = hco3.accu,stderr)
t.test(hco3.accu[-7,]$dhco3~hco3.accu[-7,]$treat)








# calculate rates of pH-compensation
data.rate<-data.frame(fishID = unique(data.comp$fishID),rate = NA, treatment = NA, SA = NA)

for (k in 1:12){
  
  # Assign treatment
  data.rate$treatment[k]<-ifelse(test = which(unique(data.rate$fishID)==data.comp$fishID[k])<7,yes = "H", no = "L")
  # Calculate rate
  temp.df <- data.comp[data.comp$fishID==data.rate$fishID[k],]
  data.rate$rate[k] <- ifelse(test = which(unique(data.rate$fishID)==data.comp$fishID[k])==9,no = summary(lm(compensation~time,temp.df))$coefficients[2,1], yes = NA)
  
  # Assign gill morphology
  data.gill<-read_excel(path = "./data_panga_gill.xlsx", sheet = "Gill") # Import data
  data.gill<-as.data.frame(data.gill)
  
  data.rate$SA[k]<-data.gill[which(data.gill$Fish==data.rate$fishID[k]),]$SA.m
  }
summary(lm(rate~SA,data.rate))


# Test if groups have different pH-compensation rates
t.test(data.rate[-9,]$rate~data.rate[-9,]$treatment)
# They dont. 

## PLOT DATA ##
# parameters to be plotted
para<-c("ph","pco2","hco3","osm","na","cl")

# Corresponding axis labels
label                    <- c(expression(bold("pH"[a])),
                              expression(bold("P"[a]*"CO"["2"]*" (mmHg)")),
                              expression(bold("[HCO"[3]^-{}*"]"[p]*" (mmol l"^"-1"*")")),
                              expression(bold("Osmolality (mOsm kg"^"-1"*")")),
                              expression(bold("[Na"^"+"*"]"[p]*" (mmol l"^"-1"*")")),
                              expression(bold("[Cl"^"-"*"]"[p]*" (mmol l"^"-1"*")")))
                              
# Generate a function that runs stats and plots each parameter as a function of salinity
fun<-function (i) {
  
  
  df <- data
  # Calculate means and standard error
  
  df$tt<-paste(ifelse(nchar(as.character(df$time))!=2,paste0("0",as.character(df$time)),as.character(df$time)),df$treatment,sep = "")
  mod.full<-lmer(paste(as.name(para[i]),"~tt+(1|fishID)"),df)
  mod.2<-lmer(paste(as.name(para[i]),"~time+(1|fishID)"),df)
  mod.1<-lmer(paste(as.name(para[i]),"~1+(1|fishID)"),df)
  summary(mod.2)
  anova(mod.2,mod.1)
  summary(mod.1)
  fit<-summary(glht(mod.full, linfct = mcp(tt = "Tukey")), test = adjusted("holm"))
  fit$test
  p<-fit$test$pvalues
  names(p)<-sub(" - ", "-", names(p))
  p<-multcompLetters(p)
  p
  pairwise<-data.frame(
    coef = fit$test$coefficients,
    sigma  = fit$test$sigma,
    t = fit$test$tstat,
    p = fit$test$pvalues
  )
  
  
  temp.p<-p
  
  if(all(temp.p$Letters=="a")==T){
    p <- rep("",10)
  }
  
  if(all(temp.p$Letters=="a")==F){
    p <- unname(p$Letters)[c(10,1:9)]
  }
  p
  
  
  a = cbind(aggregate(as.formula(paste(as.name(para[i]),"~treatment+time")),df,mean),                         # Data frame with mean and standard error over salinities
            aggregate(as.formula(paste(as.name(para[i]),"~treatment+time")),df,stderr)[,3],
            aggregate(as.formula(paste(as.name(para[i]),"~treatment+time")),df,length)[,3])
  colnames(a)<-c("Treatment","Time","Mean","stderr","n")
  a$Treatment<-as.factor(a$Treatment)
  a

  
  write.table(x = a, file = paste("./stats/",para[i],"_means.txt",sep=""))
  write.table(x = pairwise, file = paste("./stats/",para[i],"_pairwise.txt",sep=""))

  plot<-
    ggplot() +
    
    # Add individual data points and connecting lines
    lapply(
      X = 1:12, 
      FUN = function(k){
        temp.df<-data[data$fishID==unique(data$fishID)[k],];temp.df
        if(unique(temp.df$treatment)=="H"){col<-cols[1]}
        if(unique(temp.df$treatment)=="L"){col<-cols[2]}
        temp.df<-
          data.frame(
            temp.df[para[i]],
            temp.df["time"])
        temp.df<-na.omit(temp.df)
        colnames(temp.df)<-c("value","time")
        if(dim(temp.df)[1]!=0){
          geom_line(data = temp.df, mapping = aes(x = time, y = value,group = 1), color =col, alpha = 0.3)  
        }
      }) +
      
      # Add means standard error and connecting lines
      geom_point(data=a, aes(x=as.numeric(as.character(Time)), y=Mean, col = Treatment), size = 2)+
      geom_path(data=a, aes(x=as.numeric(as.character(Time)), y=Mean, col = Treatment))+
      geom_errorbar(data = a, aes(x=as.numeric(as.character(Time)), y = Mean, ymin=Mean-stderr, ymax=Mean+stderr, col = Treatment),width=0.75) +
      scale_color_manual(values = cols[1:2], labels = c("Hyperoxia", "Hypoxia"), breaks = c("H","L"))+
      
      # add pairwise comparisons if blue points are higher than red points
      lapply(
        X = 1:5, 
        FUN = function(k){
          times <- c("0","4","8","18","24")
          ps<-c(2*(k-1)+1,2*(k-1)+2)
          if(diff(a$Mean[a$Time==times[k]])>0){
            geom_text(
              mapping = aes(
                label = p[ps],
                y = c(as.numeric(a$Mean-a$stderr)[ps[1]],as.numeric(a$Mean+a$stderr)[ps[2]]),
                x = as.numeric(times[k])),
              vjust = c(1.5,-0.5), 
              color = "black" , 
              size = 2.8)
          }
          }
      ) +
      
      
      # add pairwise comparisons if red points are higher than blue points
      lapply(
        X = 1:5, 
        FUN = function(k){
          times <- c("0","4","8","18","24")
          ps<-c(2*(k-1)+1,2*(k-1)+2)
          
          # if blue points are lower than red points
          if(diff(a$Mean[a$Time==times[k]])<0){
            geom_text(
              mapping = aes(
                label = p[ps],
                y = c(as.numeric(a$Mean+a$stderr)[ps[1]],as.numeric(a$Mean-a$stderr)[ps[2]]),
                x = as.numeric(times[k])),
              vjust = c(-.5,1.5), 
              color = "black" , 
              size = 2.8)
          }
          }
        ) +
      
    # Set up plotting theme
    theme_classic()+
    theme(
      legend.position = "none",
      panel.grid.major   = element_line(colour = "grey90",size = 0.2),
      axis.text = element_text(size=6),
      axis.title = element_text(size=6, margin = margin(t=0,r = 0,b = 0,l = 10)),
      strip.background = element_blank())+
    ylab(label[i])+
    xlab(expression(bold("Time (h)")))
plot
p
temp.p
  return(plot)
}

i=1
repeat {
  fun(i)
  ggsave(filename = paste("./Figures/",para[i],".pdf",sep = ""),width = 3.5,height = 3.5)
  ggsave(filename = paste("./Figures/",para[i],".jpeg",sep = ""),width = 3.5,height = 3.5)  
  i = i+1
  if (i == length(para)+1){
    break
  }
}


pdf("./Figures/Fig. 1 - acidbase.pdf",width = 3.5,height = 10,useDingbats = F)
plot_grid(fun(1),fun(2),fun(3),nrow = 3,ncol = 1,align = 'h',labels = "AUTO",label_size = 12)
dev.off()

jpeg(filename = "./Figures/Fig. 3 - acidbase.jpeg",width = 3.5,height = 10,units = "in",res = 600)
plot_grid(fun(1),fun(2),fun(3),nrow = 3,ncol = 1,align = 'h',labels = "AUTO",label_size = 12)
dev.off()



## Davenport
dav <- function(k){
  if(k == T){
    
  
df <- data


a<-cbind(aggregate(ph~time+treatment, df, mean),
         aggregate(ph~time+treatment, df, stderr)[,3],
         aggregate(hco3~time+treatment, df, mean)[,3],
         aggregate(hco3~time+treatment, df, stderr)[,3])

colnames(a)<-c("Time","Treatment","pH","pH.Err","HCO3","HCO3.Err")


# data frame for extracellular non bicarbonate buffering effect from Damsgaard et al 2015 JEB
nbb<-data.frame(
  ph = c(7.75,7.3),
  hco3 = c(18 ,18 -18.3*(7.3-7.75))
)

cols
# Plot and save Figure 2 as PDF
ggplot() + 
  geom_line(mapping = aes(nbb$ph, y = nbb$hco3),linetype ="dashed") +
#  geom_polygon(mapping = aes(nbb2$ph, y = nbb2$hco3),alpha = 0.5) +
  #scale_x_continuous(limits = c(7.3,7.65)) +
  scale_y_continuous(limits = c(15,45)) +
  geom_point(a, mapping = aes(x = pH, y = HCO3,col = Treatment)) +
  geom_path(a, mapping = aes(x = pH, y = HCO3, col = Treatment)) +
  geom_errorbar(width = 0.015,aes(x = a[,3], ymin=a[,5]-a[,6], ymax=a[,5]+a[,6], col = a$Treatment)) +
  geom_errorbarh(height = 0.6, aes(y = a[,5], xmin=a[,3]-a[,4], xmax=a[,3]+a[,4], col = a$Treatment)) +
  scale_color_manual(values = cols[1:2], breaks = c("H", "L"), labels = c("Hyperoxia", "Hypoxia"))+
  theme_classic()+
  theme(legend.position=c(0.2,0.95),
        panel.grid.major   = element_line(colour = "grey90",size = 0.2),
        legend.text = element_text(size=6, face = "bold"),
        legend.title = element_blank(),
        axis.text = element_text(size=6),legend.background = element_rect(fill = "transparent"),
        axis.line = element_line(size = 0.5),
        axis.title = element_text(size=6, margin = margin(t=0,r = 0,b = 0,l = 10)),
        strip.background = element_blank(),
        plot.margin = unit(c(0.1, 0, 0, 0.1), "cm"))+
  geom_text(mapping = aes(x = a[,3], y = a[,5], label=c("0 h","4 h","8 h","18 h","24 h","0 h","4 h","8 h","18 h","24 h")), 
            hjust=-.5, vjust = 1.5,
            color="black",
            #position = position_dodge(0.9), 
            size=2.8)+
  annotate("text",size=2.8, x = c(7.74,7.8,7.75,7.6), y = c(15,33,45,45), vjust = 0, hjust = 0.5, label = c("10","20","30","40"),fontface = "bold")+
  stat_function(fun = function(pH) 10*10^(pH-(4.83+0.17*pH))*aco2, linetype = "dotted", colour="black")+
  stat_function(fun = function(pH) 20*10^(pH-(4.83+0.17*pH))*aco2, linetype = "dotted", colour="black")+
  stat_function(fun = function(pH) 30*10^(pH-(4.83+0.17*pH))*aco2, linetype = "dotted", colour="black")+
  stat_function(fun = function(pH) 40*10^(pH-(4.83+0.17*pH))*aco2, linetype = "dotted", colour="black")+
  labs(x = expression(bold("pH"[a])),
       y = expression(bold("[HCO"[3]^-{}*"]"[p]*" (mmol l"^"-1"*")"))) -> plot
return(plot)
  }
  }
dav(k=T)


### Rate of pH compensation
df<-data.comp

a = cbind(aggregate(compensation~treatment+time,df,mean),                         # Data frame with mean and standard error over salinities
          aggregate(compensation~treatment+time,df,stderr)[,3],
          aggregate(compensation~treatment+time,df,length)[,3])
colnames(a)<-c("Treatment","Time","Mean","stderr","n")
a$Treatment<-as.factor(a$Treatment)



rate.panel1<-
  ggplot() +
  
    # Add means standard error and connecting lines
    
    
    # Add individual data points and connecting lines
    lapply(
      X = 1:12, 
      FUN = function(k){
        temp.df<-data.comp[data.comp$fishID==unique(data.comp$fishID)[k],];temp.df
        if(unique(temp.df$treatment)=="H"){col<-cols[1]}
        if(unique(temp.df$treatment)=="L"){col<-cols[2]}
        temp.df<-
          data.frame(
            temp.df["compensation"],
            temp.df["time"])
        temp.df<-na.omit(temp.df)
        colnames(temp.df)<-c("value","time")
        if(dim(temp.df)[1]!=0){
          geom_line(data = temp.df, mapping = aes(x = time, y = value,group = 1), color =col, alpha = 0.3)  
        }
      })+
    geom_point(data=a, aes(x=Time, y=Mean, col = Treatment))+
    geom_path(data=a, aes(x=Time, y=Mean, col = Treatment))+
    geom_errorbar(data = a, aes(x = Time, y = Mean, ymin=Mean-stderr, ymax=Mean+stderr, col = Treatment),width=0.75) +
    scale_color_manual(values = cols[1:2])+
    # Set up plotting theme
    theme_classic()+
    theme(
      legend.position = "none",
      panel.grid.major   = element_line(colour = "grey90",size = 0.2),
      axis.text = element_text(size=6),
      axis.title = element_text(size=6, margin = margin(t=0,r = 0,b = 0,l = 10)),
      strip.background = element_blank())+
    ylab(expression(bold("Percent pH-compensation")))+
    xlab(expression(bold("Time (h)")))
rate.panel1
rate.panel2<-
  ggplot()+
    geom_point(mapping = aes(x = SA, y = 100*rate, col = treatment), data = data.rate)+
    scale_color_manual(values = cols[1:2])+
    # Set up plotting theme
    theme_classic()+
    theme(
      legend.position = "none",
      panel.grid.major   = element_line(colour = "grey90",size = 0.2),
      axis.text = element_text(size=6),
      axis.title = element_text(size=6, margin = margin(t=0,r = 0,b = 0,l = 10)),
      strip.background = element_blank())+
    ylab(expression(bold("Rate of pH-compensation (% hr"^"-1"*")")))+
    xlab(expression(bold("Gill surface area (mm"^"2"*" kg"^"-1"*")")))

pdf("./Figures/Fig. 4 - compensation.pdf",width = 3.5,height = 7,useDingbats = F)
plot_grid(rate.panel1,rate.panel2,nrow = 2,ncol = 1,align = 'v',labels = "AUTO",label_size = 12)
dev.off()

jpeg("./Figures/Fig. 4 - compensation.jpeg",width = 3.5,height = 7,units = "in",res = 600)
plot_grid(rate.panel1,rate.panel2,nrow = 2,ncol = 1,align = 'v',labels = "AUTO",label_size = 12)
dev.off()


  plot_grid(
    plot_grid(fun(1),fun(2),fun(3),nrow = 1,ncol = 3,align = 'h',labels = c("A","B","C"),label_size = 12),
    plot_grid(dav(k=T),rate.panel1,rate.panel2,nrow = 1,ncol = 3,align = 'h',labels = c("D","E","F"),label_size = 12,rel_widths = c(2,1,1)),
    ncol = 1, nrow = 2)

  
  
  plot_grid(
    plot_grid(dav(k=T),labels = "A", label_size = 12),
    plot_grid(fun(1),fun(2),fun(3),rate.panel1,nrow = 2,ncol = 2,align = 'hv',labels = c("B","C","D","E"),label_size = 12),
    ncol = 2, nrow = 1)

jpeg("./Figures/Figure 2.jpeg",width = 12,height = 18,units = "cm",res = 600)
  plot_grid(
  plot_grid(dav(k=T),labels = "A", label_size = 12),
  plot_grid(fun(1),fun(2),fun(3),rate.panel1,nrow = 2,ncol = 2,align = 'hv',labels = c("B","C","D","E"),label_size = 12),ncol = 1, nrow = 2)
dev.off()
## Gill morphology
colnames(data.gill)

para.gill <- c("SA.m","HM","ADF")
label.gill<-c(expression(bold("Gill surface area (mm"^"2"*" kg"^"-1"*")")),
          expression(bold("Harmonic mean (µm)")),
          expression(bold("ADF (cm"^"2"*" µm"^"-1"*" g"^"-1"*")")))


# Generate a function that runs stats and plots each parameter as a function of salinity
fun.gill<-function (i) {
  
  # Calculate means and standard error
  df<-data.gill
  
  df$treatment <- as.factor(df$Treatment)
  
  a = cbind(aggregate(as.formula(paste(as.name(para.gill[i]),"~treatment")),df,mean),                         # Data frame with mean and standard error over salinities
            aggregate(as.formula(paste(as.name(para.gill[i]),"~treatment")),df,stderr)[,2],
            aggregate(as.formula(paste(as.name(para.gill[i]),"~treatment")),df,length)[,2])
  colnames(a)<-c("Treatment","Mean","stderr","n")
  a$Treatment<-as.factor(a$Treatment)
  write.table(x = a, file = paste("./stats/",para.gill[i],"_means.txt",sep=""))
  
  ttest<-t.test(as.formula(paste(as.name(para.gill[i]),"~treatment")),data = df)
  write.table(data.frame(t = ttest$statistic, df = ttest$parameter, ttest$p.value),file = paste("./stats/",para.gill[i],"_ttest.txt",sep=""))
  
  testlabel<-rep(NA,2)
  if(ttest$p.value<0.05){
    testlabel <- c("a","b")
  }
  
  raw<-data.frame(df[para.gill[i]], df["Treatment"])
  colnames(raw)<-c("value","treatment")
  
  plot<-
    ggplot() +
    scale_x_discrete(breaks = c("H","L"), labels = c("Hyperoxia","Hypoxia"))+
    scale_color_manual(breaks = c("H","L"), values = cols[1:2])+
    
    geom_point(data=a, aes(x = Treatment, y=Mean), size = 2)+
    geom_errorbar(data = a, aes(x = Treatment, y = Mean, ymin=Mean-stderr, ymax=Mean+stderr), width=0.25) +
    geom_point(mapping = aes(x = treatment, y = value, col = treatment),data = raw) +
    geom_text(
      mapping = aes(
        label = testlabel,
        y = c(max(raw[raw$treatment=="H",]$value,na.rm = T),max(raw[raw$treatment=="L",]$value,na.rm = T)),
        x = c("H","L")),
      vjust = c(-1), 
      color = "black" , 
      size = 2.8)+
    ylim(c(0.8*min(raw$value,na.rm = T),1.2*max(raw$value,na.rm = T))) +
  
    # Set up plotting theme
    theme_classic()+
    theme(
      legend.position = "none",
      panel.grid.major   = element_line(colour = "grey90",size = 0.2),
      axis.text = element_text(size=6),
      axis.title = element_text(size=6, margin = margin(t=0,r = 0,b = 0,l = 10)),
      strip.background = element_blank())+
    ylab(label.gill[i])+
    xlab("")
  plot
  
  return(plot)
}

i=1
repeat {
  fun.gill(i)
  ggsave(filename = paste("./Figures/",para.gill[i],".pdf",sep = ""),width = 2,height = 3.5)
  ggsave(filename = paste("./Figures/",para.gill[i],".jpeg",sep = ""),width = 2,height = 3.5)  
  i = i+1
  if (i == length(label.gill)+1){
    break
  }
}



plot_grid(fun.gill(1),fun.gill(2),fun.gill(3),nrow = 3,ncol = 1,align = 'h',labels = "AUTO",label_size = 12,hjust = -2)
pdf("./Figures/Fig. 1 - gill morphology.pdf",width = 2.75,height = 5,useDingbats = F)
plot_grid(fun.gill(1),fun.gill(2),fun.gill(3),nrow = 3,ncol = 1,align = 'v',labels = c("C","D","E"),label_size = 12, hjust = -4.5)
dev.off()

jpeg("./Figures/Fig. 1 - gill morphology.jpeg",width = 3.5,height = 10,units = "in",res = 600)
plot_grid(fun.gill(1),fun.gill(2),fun.gill(3),nrow = 3,ncol = 1,align = 'h',labels = "AUTO",label_size = 12)
dev.off()
