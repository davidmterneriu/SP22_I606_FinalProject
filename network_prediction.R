################################################################################
# BlueWater Compair Network Prediction
################################################################################
rm(list=ls())

library(tidyverse)
library(ggplot2)
library(readr)
library(readxl)
library(fixest)
library(gganimate)
library(gifski)
library(transformr)


bw_edge<-read_csv("~/Desktop/IU/INFO/BW_world/bw_cap_edge_df.csv")
bw_ports<-read_csv("~/Desktop/IU/INFO/BW_world/port_flows.csv")

bw_port_df=rbind(bw_ports%>%select(time_m,port=to,cap=flow),
      bw_ports%>%select(time_m,port=from,cap=flow))%>%
  group_by(time_m,port)%>%
  summarise(flow=sum(cap))%>%
  ungroup()%>%
  mutate(flow_ln=log10(flow+1))%>%
  mutate(year=floor(time_m),
         month=(time_m-year)*12+1,
         mon.abb=month.abb[month],
         time_str=paste0(year," ",mon.abb))

annual_sum=read_csv("~/Desktop/IU/INFO/BW_world/annual_sum.csv")

annual_sum=annual_sum%>%
  mutate(mon.abb=month.abb[month_n],
         time_str=paste0(year," ",mon.abb))

g1=annual_sum%>%
  ggplot(aes(x=time_m,y=tot_cap/(10^9),group=1))+
  geom_line()+
  geom_point()+
  theme_bw()+
  labs(subtitle = "Time: {frame_along}",y="TEU Capacity (billions)",x="Time")+
  transition_reveal(time_m,keep_last=FALSE)+
  ease_aes("linear")+
  shadow_mark()
animate(g1, width = 500, height = 500, fps = 2, duration = 20,
        end_pause = 10, res = 100,renderer = gifski_renderer())
anim_save("~/Desktop/IU/INFO/BW_world/tot_teu.gif")




#install.packages("magick")


port_country_cords <- read_csv("~/Desktop/IU/INFO/BW_world/port_country_cords.csv")

bw_port_df=bw_port_df%>%inner_join(port_country_cords)

graph2=ggplot(data=bw_port_df,aes(x=lon,y=lat,size=flow))+
  geom_point()+
  labs(size="Port TEU")+
  ggthemes::theme_map()



graph1.animation = graph2 +
  transition_time(time_m) +
  labs(subtitle = "Year: {frame_time}")+
  ease_aes('linear')


animate(graph1.animation, width = 500, height = 500, fps = 2, duration = 20,
        end_pause = 10, res = 100,renderer = gifski_renderer())
anim_save("~/Desktop/IU/INFO/BW_world/port_teu.gif")



library(magick)

a_mgif <- image_read("~/Desktop/IU/INFO/BW_world/port_teu.gif")
b_mgif <- image_read("~/Desktop/IU/INFO/BW_world/tot_teu.gif")

new_gif <- image_append(c(a_mgif[1], b_mgif[1]))
for(i in 2:30){
  combined <- image_append(c(a_mgif[i], b_mgif[i]))
  new_gif <- c(new_gif, combined)
}
new_gif
anim_save("~/Desktop/IU/INFO/BW_world/combined_port_map.gif")


dist_cepii <- read_excel("~/Desktop/IU/INFO/BW_world/dist_cepii.xls")
dist_cepii=dist_cepii[,-c(11:14)]

port_grid=expand.grid(port_i=port_country_cords$port,port_j=port_country_cords$port)%>%
  as.data.frame()%>%
  inner_join(port_country_cords,by=c("port_i"="port"))%>%
  inner_join(port_country_cords,by=c("port_j"="port"))%>%
  filter(port_i!=port_j)

port_grid=port_grid%>%
  mutate(dist=geosphere::distHaversine(p1=cbind(lon.x,lat.x),p2=cbind(lon.y,lat.y),
                                       r=3958.8))


port_grid=port_grid%>%
  select(port_i,port_j,iso3c_i=country.x,iso3c_j=country.y,dist)%>%
  left_join(dist_cepii,by=c("iso3c_i"="iso_o","iso3c_j"="iso_d"))

port_grid[is.na(port_grid)]<-0

bw_ports1=bw_ports%>%
  arrange(time_m)%>%
  group_by(to,from)%>%
  mutate(lead_flow=lead(flow,1))%>%
  ungroup()%>%
  filter(is.na(lead_flow)==F)

################################################################################
# Creating Test Data
################################################################################

max_t=bw_ports1$time_m%>%max()


bw_ports1=bw_ports1%>%inner_join(port_grid,by=c("to"="port_i","from"="port_j"))

bw_ports1$year=floor(bw_ports1$time_m)
bw_ports1_train=filter(bw_ports1,time_m<max_t)

bw_edge_train=filter(bw_edge,time_id<=2021.70)


bw_ports1_test=filter(bw_ports1,time_m==max_t)
bw_edge_test=filter(bw_edge,time_id>2021.70)%>%
  filter(tot_vol_lead>0)



bw_ports1_test%>%colnames()

grav_test=feols(data=filter(bw_ports1_train,lead_flow>0,dist>0),
                log(lead_flow)~log(dist)+shared_carrier+contig+col45+comlang_off+smctry|
                  to^year+from^year,combine.quick=FALSE)

grav_fe=fixest::fixef(grav_test)
to_year_df=data.frame(theta_it=grav_fe$`to^year`)
to_year_df$to_year<-rownames(to_year_df)
rownames(to_year_df)<-c()

from_year_df=data.frame(eta_it=grav_fe$`from^year`)
from_year_df$from_year<-rownames(from_year_df)
rownames(from_year_df)<-c()





dist_imped=function(alpha){
  #browser()
  temp_df=filter(bw_ports1_train,lead_flow>0,dist>0)
  
  temp_df=temp_df%>%
    mutate(jc_0_a=(jc+1)/dist^alpha,
           jc_1_a=(jc1+1)/dist^alpha,
           TEU_0=TEU_i+1,
           TEU_1=TEU_j+1)%>%
    group_by(time_m)%>%
    mutate(r1_0_a=log(10^(-6)+jc_0_a/sum(jc_0_a,na.rm=T)),
           r1_1_a=log(10^(-6)+jc_1_a/sum(jc_1_a,na.rm=T)))%>%
    ungroup()
  
  net_mod=lm(data=temp_df,
             log(lead_flow)~log(TEU_0)+log(TEU_1)+log(degree_i)+log(degree_j)+log(close_i)+log(close_j)+
                r1_0_a+r1_1_a+between_i+between_j+coreness_i+coreness_j+ecc_i+ecc_j)
  error=mean(resid(net_mod)^2)
  return(error)
}
dist_imped(100)

opt_alpha=optim(par=1,fn=dist_imped,method = "Brent",
      lower=c(0),upper = c(100),hessian = T)


alpha_1=opt_alpha$par

bw_ports1_train2=filter(bw_ports1_train,lead_flow>0,dist>0)%>%
  mutate(jc_0_a=(jc+1)/dist^alpha_1,
       jc_1_a=(jc1+1)/dist^alpha_1,
       TEU_0=TEU_i+1,
       TEU_1=TEU_j+1)%>%
  group_by(time_m)%>%
  mutate(r1_0_a=log(10^(-6)+jc_0_a/sum(jc_0_a)),
         r1_1_a=log(10^(-6)+jc_1_a/sum(jc_1_a)))%>%
  ungroup()


bw_ports1_test2=bw_ports1_test%>%mutate(jc_0_a=(jc+1)/dist^alpha_1,
                       jc_1_a=(jc1+1)/dist^alpha_1,
                       TEU_0=TEU_i+1,
                       TEU_1=TEU_j+1)%>%
  group_by(time_m)%>%
  mutate(r1_0_a=log(10^(-6)+jc_0_a/sum(jc_0_a)),
         r1_1_a=log(10^(-6)+jc_1_a/sum(jc_1_a)))%>%
  ungroup()%>%
  filter(lead_flow>0)

bw_ports1_test2=bw_ports1_test2%>%
  mutate(to_year=paste0(to,"_",year),
         from_year=paste0(from,"_",year))%>%
  inner_join(to_year_df)%>%
  inner_join(from_year_df)






net_test=lm(data=bw_ports1_train2,
               log(lead_flow)~log(TEU_0)+log(TEU_1)+log(degree_i)+log(degree_j)+log(close_i)+log(close_j)+
              r1_0_a+r1_1_a+between_i+between_j+coreness_i+coreness_j+ecc_i+ecc_j)

library(shapley)

reg <- function(factors, dv, data) {
  if (length(factors) == 0) return(0)
  formula <- paste(dv, "~", paste(factors, collapse = "+"))
  m <- lm(formula, data = data)
  summary(m)$r.squared
}



grav_test=feols(data=filter(bw_ports1_train,lead_flow>0,dist>0),
                log(lead_flow)~log(dist)+shared_carrier+contig+col45+comlang_off+smctry|
                  to^year+from^year,combine.quick=FALSE)

reg2 <- function(factors,  data) {
  if (length(factors) == 0) return(0)
  form <- paste("log(lead_flow)", "~", paste(factors, collapse = "+"),"|to^year+from^year")
  m <- feols(data = data,as.formula(form))
  res=cor(fitted(m),log(data$lead_flow))^2
  return(res)
}

reg2(data=filter(bw_ports1_train,lead_flow>0,dist>0),
     factors =names(grav_test$coefficients) )


pred_names=colnames(net_test$model)
pred_names=pred_names[-1]

if(FALSE){
shape_res=shapley::shapley(reg, pred_names ,silent = TRUE, dv="log(lead_flow)", data = bw_ports1_train2)

shape_res%>%
  ggplot(aes(x=reorder(factor,-value),y=value))+
  geom_bar(stat="identity",color="black")+
  geom_text(aes(label=round(value,3)),vjust=-0.2)+
  labs(x="Predictor",y="Contribution to R2",
       title = "Network Features Model: Shapley Value Decomposition")+
  theme_bw()

sum(shape_res$value)


shape_res2=shapley::shapley(reg2, names(grav_test$coefficients) ,silent = FALSE, 
                            data = filter(bw_ports1_train,lead_flow>0,dist>0))

sum(shape_res2$value)



}

edge_mod=lm(data=bw_edge_train%>%select(-c(source,target,time_id))%>%filter(tot_vol_lead>0),
            log(tot_vol_lead)~.)


library(caret)
library(glmnet)
set.seed(42)
cv_5 = trainControl(method = "cv", number = 50)

net_feat_elnet = train(
  log(lead_flow)~log(degree_i)+log(degree_j)+log(close_i)+log(close_j)+
    r1_0_a+r1_1_a+between_i+between_j+coreness_i+coreness_j+ecc_i+ecc_j,
  data=bw_ports1_train2,
  method = "glmnet",
  preProcess = c("center", "scale"),
  trControl = cv_5
)



net_elnet = train(
  log(tot_vol_lead)~., 
  data=bw_edge_train%>%select(-c(source,target,time_id))%>%filter(tot_vol_lead>0),
  method = "glmnet",
  preProcess = c("center", "scale"),
  trControl = cv_5
)


get_best_result = function(caret_fit) {
  best = which(rownames(caret_fit$results) == rownames(caret_fit$bestTune))
  best_result = caret_fit$results[best, ]
  rownames(best_result) = NULL
  best_result
}

rbind(get_best_result(net_feat_elnet),
      get_best_result(net_elnet))%>%
  kableExtra::kable(format = "latex",digits = 3,booktabs=T)

net_elnet$modelInfo$predict()

net_elnet$bestTune


texreg::texreg(grav_test,digits = 3,stars = c(0.01,0.05,0.1))

################################################################################
#Fit
################################################################################

cor(predict(grav_test,newdata = bw_ports1_test2),log(bw_ports1_test2$lead_flow))^2


cor(predict(net_feat_elnet,bw_ports1_test2),log(bw_ports1_test2$lead_flow))^2
cor(predict(net_elnet, bw_edge_test),log(bw_edge_test$tot_vol_lead))^2


compare_df=rbind(data.frame(pred=predict(grav_test,newdata =filter(bw_ports1_train,lead_flow>0,dist>0)),
           obs=filter(bw_ports1_train,lead_flow>0,dist>0)%>%select(obs=lead_flow),
           model="Gravity",
           type="Train"),
      data.frame(pred=predict(grav_test,newdata = bw_ports1_test2)%>%as.numeric(),
                 obs=bw_ports1_test2%>%select(obs=lead_flow),
                 model="Gravity",
                 type="Test"),
      #Net Features
      data.frame(pred=predict(net_feat_elnet,newdata =filter(bw_ports1_train2,lead_flow>0,dist>0)),
                 obs=filter(bw_ports1_train,lead_flow>0,dist>0)%>%select(obs=lead_flow),
                 model="Network:Features",
                 type="Train"),
      data.frame(pred=predict(net_feat_elnet,newdata = bw_ports1_test2),
                 obs=bw_ports1_test2%>%select(obs=lead_flow),
                 model="Network:Features",
                 type="Test"),
      #Embedd
      data.frame(pred=predict(net_elnet,newdata =bw_edge_train),
                 obs=filter(bw_edge_train)%>%select(obs=tot_vol_lead),
                 model="Network:Embed",
                 type="Train"),
      data.frame(pred=predict(net_elnet,newdata=bw_edge_test),
                 obs=bw_edge_test$tot_vol_lead,
                 model="Network:Embed",
                 type="Test"))

compare_df%>%
  filter(obs>0)%>%
  mutate(obs=log(obs))%>%
  group_by(model,type)%>%
  summarise(R2=cor(pred,obs)^2)%>%
  ggplot(aes(x=model,y=R2,fill=model))+
  geom_bar(stat="identity",position = "dodge")+
  facet_wrap(~type)+
  labs(x="",y="R2",fill="Model:")+
  geom_text(aes(label=round(R2,2)),vjust=-0.2)+
  ggthemes::theme_clean()+
  scale_y_continuous(n.breaks = 10)+
  theme(legend.position = "top",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggthemes::scale_fill_calc()


compare_df%>%
  filter(obs>0)%>%
  mutate(obs=log(obs))%>%
  group_by(model,type)%>%
  summarise(R2=mean((pred-obs)^2))%>%
  ggplot(aes(x=model,y=R2,fill=model))+
  geom_bar(stat="identity",position = "dodge")+
  facet_wrap(~type)+
  labs(x="",y="RMSE",fill="Model:")+
  geom_text(aes(label=round(R2,2)),vjust=-0.2)+
  ggthemes::theme_clean()+
  theme(legend.position = "top",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggthemes::scale_fill_calc()


compare_df%>%
  filter(obs>0)%>%
  mutate(obs=log(obs),
         pe_error=(obs-pred))%>%
  ggplot(aes(x=model,y=pe_error,color=model))+
  geom_boxplot()+
  facet_wrap(~type)



