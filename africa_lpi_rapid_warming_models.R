library(broom)
library(coefplot)
library(dplyr)
library(ggplot2)
library(lme4)
library(MuMIn)

source("rsquaredglmm.R")

df <- read.csv("africa_lpi_rapid_warming_data.csv")
#df <- read.csv("africa_lpi_rapid_warming_data_norsq.csv")
df$bodymass_g <- 10^df$Log_Body_Mass_g
#Identifying sites with more than one population trend

sp_dups <- df %>% 
  group_by(Longitude, Latitude) %>% 
  summarise(n = n()) 

sp_dups$loc_id <- 1:length(sp_dups$Longitude)

sp_dups_df<-merge(sp_dups, df, by=c("Longitude","Latitude"))

parm_df <- sp_dups_df %>% 
  dplyr::select("ID","mean_annual_temp_change", "annual_crop_gras_change", "bodymass_g")

parm_scale <- scale(parm_df[,c("mean_annual_temp_change", "annual_crop_gras_change", "bodymass_g")])       #use the scaling factors at the bottom of these to scale the rasters

parm_id <- parm_df[,"ID"]

parm_df_scale <- data.frame(parm_id, parm_scale)

colnames(parm_df_scale) <- c("ID","mean_temp_scale", "luc_change_scale", "bodymass_scale")

dt <- merge(sp_dups_df, parm_df_scale, by="ID")

nrow(dt)


m0 <- lmer(lambda_mean ~ luc_change_scale + mean_temp_scale + luc_change_scale:mean_temp_scale+bodymass_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m0a <- lmer(lambda_mean ~ luc_change_scale + mean_temp_scale + bodymass_scale+(1|Binomial) + (1|loc_id),data=dt, REML=F)

m0b <- lmer(lambda_mean ~ luc_change_scale+ bodymass_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m0c <- lmer(lambda_mean ~ mean_temp_scale + bodymass_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m0d <- lmer(lambda_mean ~ bodymass_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m1 <- lmer(lambda_mean ~ luc_change_scale + mean_temp_scale + luc_change_scale:mean_temp_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m1a <- lmer(lambda_mean ~ luc_change_scale + mean_temp_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m1b <- lmer(lambda_mean ~ luc_change_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

m1c <- lmer(lambda_mean ~ mean_temp_scale + (1|Binomial) + (1|loc_id), data=dt, REML=F)

mnull <-lmer(lambda_mean ~ 1 + (1|Binomial) + (1|loc_id), data=dt, REML=F)


###Model Comparison

msAICc <- model.sel(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)

msAICc$model<-rownames(msAICc)
msAICc<-data.frame(msAICc)
msAICc

#Rsq
models_list<-list(m0,m0a,m0b,m0c,m0d,m1,m1a,m1b,m1c,mnull)
modelsR<-lapply(models_list,rsquared.glmm)
modelsRsq <- matrix(unlist(modelsR), ncol=6, byrow=T)

modelsRsq

#mav<-model.avg(models_list, subset = delta <= 2)
mav<-model.avg(models_list, subset = cumsum(weight) <= .95)

smav<-summary(mav)

coef_av<-smav$coefmat.full[,"Estimate"]
coef_df<-data.frame(coef_av)
coef_df$lowCI<-confint(mav)[,1]
coef_df$highCI<-confint(mav)[,2]
coef_df

coef_pcnt<-data.frame(((10^coef_df) - 1)*100) #converting to annual percentage change
coef_pcnt$Var_name <- row.names(coef_pcnt)
coef_pcnt$Class <- "Mammalia"

round(coef_pcnt[,1:3],3)


p1 <- ggplot(coef_pcnt) +
  geom_linerange(
    aes(x = Var_name, ymin = lowCI, ymax = highCI),
    lwd = 2.5,
    position = position_dodge(width = 2 / 3)
  ) +
  geom_pointrange(
    aes(
      x = Var_name,
      y = coef_av,
      ymin = lowCI,
      ymax = highCI
    ),
    lwd = 2,
    position = position_dodge(width = 2 / 3),
    shape = 21,
    fill = "White"
  ) +
  labs(y = "Average Annual Rate of\nPopulation Change (%)", x = "") +
  theme_bw() +
  theme(
    legend.position = "none",
    text = element_text(size = 20),
    axis.text.x = element_text(size = 20),
    axis.title.x = element_text(margin = unit(c(5, 0, 0, 0), "mm")),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  scale_color_manual(values = c("black", "black")) +
  geom_hline(yintercept = 0, linetype = 2) + 
  geom_vline(xintercept = 1.5, linetype = 2) + 
  scale_x_discrete(labels = c("(Intercept)" = "Intercept",                                                        
                              "mean_temp_scale" = "Rate of climate \nwarming \n(RCW)",                     
                              "luc_change_scale" = "Rate of \nconversion to \nanthropogenic \nland use \n(RCA)", 
                              "luc_change_scale:mean_temp_scale" = "RCW:RCA",
                              "bodymass_scale" = "Body mass")
  )
  

p1

png(filename="african_mammals_lpi.png",  units="in", width=12, height=8, pointsize=12, res=1000)

print(p1)

dev.off()

###interaction plot - difficult to visualise without random effects

scale_spec <- function(x, sc, cen){
  
  sc_out <- (x-cen) / sc 
  return(sc_out)
}


scal_att<- attributes(parm_scale)

nd <- expand.grid(mean_temp = c(0.03, 0.06), luc_change = 0, bodymass = mean(dt$bodymass_g))
nd$mean_temp_scale <- scale_spec(nd$mean_temp, sc = scal_att$`scaled:scale`[1], scal_att$`scaled:center`[1])
nd$luc_change_scale <- scale_spec(nd$luc_change, sc = scal_att$`scaled:scale`[2], scal_att$`scaled:center`[2])
nd$bodymass_scale <- scale_spec(nd$bodymass, sc = scal_att$`scaled:scale`[3], scal_att$`scaled:center`[3])

lm_out <- predict(mav,nd, re.form = NA)
((10^lm_out) - 1)*100







new_d <- expand.grid(luc_change=seq(-0.02:0.03, by = 0.005),mean_temp = seq(-0.1,0.15, by = 0.025), bodymass = 808764.4)

new_d$mean_temp_scale <- scale_spec(new_d$mean_temp, sc = scal_att$`scaled:scale`[1], scal_att$`scaled:center`[1])
new_d$luc_change_scale <- scale_spec(new_d$luc_change, sc = scal_att$`scaled:scale`[2], scal_att$`scaled:center`[2])
new_d$bodymass_scale <- scale_spec(new_d$bodymass, sc = scal_att$`scaled:scale`[3], scal_att$`scaled:center`[3])

new_d$lambda_mean <- predict(mav,new_d, re.form = NA)

ggplot(dt,aes(x=mean_temp_scale,y=lambda_mean,color=luc_change_scale))+
  geom_point()+
  geom_line(data=new_d,aes(group = mean_temp_scale)) +
  scale_color_continuous(low="green",high="red")



unscale_spec <- function(x, sc, cen){
  unsc_out <- (x + cen)* sc 
  return(unsc_out)
}