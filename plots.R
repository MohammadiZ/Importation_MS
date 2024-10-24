## Winter 2022
## Project: Importation Risk
## Zahra Mohammadi
#-----------------------------
source("Epimodel.R")
## IMP and dfplot data frames will be used for plotting
## IMP : output epi model Canadian provinces, rotational workers, regular travelers, INT
## dfplot : included best model, epi model and reported CA and INT

###---------- make monthly version dataframe ------------
dfplt$year  <- strftime(dfplt$date, "%Y")
dfplt$month <- strftime(dfplt$date, "%m")

dfplt.month <- dfplt %>%
  group_by(year, month) %>%
  dplyr::summarize(CA = sum(Est.Dom), INT = sum(Est.Int),
                   bCA = sum(bestDom), bInt= sum(bestInt)
  ) %>%
  as.data.frame()

#monthly reported case
Rep  <- select(Travel_case, c(date, Tot.Int,Tot.Dom )) 
Rep$year  <- strftime(Rep$date, "%Y")   
Rep$month <- strftime(Rep$date, "%m")
Rep.month <- Rep %>%
  group_by( year, month) %>%
  dplyr::summarize( RInt= sum(Tot.Int),
                    RCA = sum(Tot.Dom)
  ) %>%
  as.data.frame()

dfplt.month  <- dfplt.month  %>% left_join(Rep.month , by ="month")
#add column date to it 
dfplt.month <- dfplt.month %>%
  mutate(date = with(., sprintf("%s-%02s", year.x, month)))
dfplt.month$date <- as.Date(paste(dfplt.month$date,"-01",sep=""))  # convert to date

##-------------------- plots ---------------------
### 1. estimated infected - bar plots
df1.sub <- select(IMP, c("report" , "RE.Int", "REW.Ca" ))
df1.sub <- df1.sub %>% rename(INT = RE.Int, CA = REW.Ca)
mdf1.sub <- melt(df1.sub, id=c("report"))

m1 <- ggplot() +
  geom_bar(data=mdf1.sub , aes(x=report, y=value, fill=variable), stat="identity") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(color=" ", 
       title="Mechanisitic model predictions",
       y="\n Reported travel-related cases (daily)",
       x="") +
  guides(fill = guide_legend(title = "Departure Origin"))+
  theme(legend.position="right") +
  scale_fill_manual(values=c(cpalete[1], cpalete[2])) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)), 
                          legend.position="bottom" ,
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))


### 2. infected percentage rw and r in Canada and International
df2.sub <- select(IMP, c("report", "RE.Int", "RE.Ca", "RW", "Total"))

for (i in 1:nrow(df2.sub)){
  for (j in 2:length(df2.sub)){
    df2.sub[i,j] <- (df2.sub[i,j]*100)/(df2.sub$Total[i])
  }
}
df2.sub <- select(df2.sub, -c("Total")) %>% rename(CA.r = RE.Ca, CA.rw = RW, INT = RE.Int)
mdf2.sub <- melt(df2.sub, id=c("report"))

m2 <- ggplot() +
  geom_bar(data=mdf2.sub , aes(x=report, y=value, fill=variable), stat="identity") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(color=" ", 
       title="Mechanistic model predictions",
       y="\n Reported travel-related cases (daily %)",
       x="") +
  guides(fill = guide_legend(title = "Type of Travelers"))+
  theme(legend.position="right") +
  scale_fill_manual(values=c(cpalete[1], cpalete[2], "#FFCC63")) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.3)), 
                          legend.position="bottom" ,
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))+ 
  scale_y_continuous(labels = function(x) paste0(x*1, "%"))


### 3. comparison models plots
### 3.1. daily International
df1.Int <- select(dfplt, c(date, bestInt, Est.Int))
df1.Int <- melt(df1.Int, id=c("date"))

plt1 <- c( "Best model"=cpalete[8], "Mechanisitic model" = cpalete[7],"Data" = cpalete[1])
p1 <- ggplot() +
  geom_line(data=dfplt, aes(x=date, y= Rep.Int, color='Data'), size=1) +
  geom_bar(data=df1.Int, aes(x=date, y=value, fill=variable), stat="identity", position=position_dodge(), alpha=.8 ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11))+
  labs(x="", y="\n Reported travel-related cases (daily)",
       title="Comparison of model predictions to reported travel-related cases \n International") +
  scale_fill_manual(name = "", labels = c( "Best model", "Mechanistic model","Data"), values = c( cpalete[8],cpalete[7], cpalete[1])) +
  scale_color_manual(values = plt1) +
  guides(color = guide_legend(title = ""), size=rel(1.4)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)),
                          legend.position="bottom",
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))


#### 3.2. daily Canada
df1.Dom <- select(dfplt, c(date, bestDom, Est.Dom))
df1.Dom <- melt(df1.Dom, id=c("date"))

plt2 <- c( "Best model"=cpalete[8],"Mechanistic model" = cpalete[7], "Data" = cpalete[2])
p2 <- ggplot() +
  geom_line( data=dfplt, aes(x=date, y= Rep.Dom, color='Data'), size=1) +
  geom_bar(data=df1.Dom, aes(x=date, y=value, fill=variable), stat="identity", position=position_dodge(), alpha=.8) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11))+
  labs(x="", y="\n Reported travel-related cases (daily)",
       title="Comparison of model predictions to reported travel-related cases \n Canada") +
  scale_fill_manual(name = "", labels = c( "Best-model", "Mechanistic-model","Reported"), values = c( cpalete[8],cpalete[7], cpalete[2])) +
  scale_color_manual(values = plt2) +
  guides(color = guide_legend(title = ""), size=rel(1.4)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)),
                          legend.position="bottom",
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))


### 3.3. monthly International
df2.Int <- select(dfplt.month, c(date, bInt, INT))
df2.Int <- melt(df2.Int, id=c("date"))

p3 <- ggplot()+
  geom_bar(data=df2.Int, aes(x=date, y=value, fill=variable), stat="identity", position=position_dodge(),alpha=.8) +
  geom_line(data= dfplt.month, aes(x=date, y=RInt),stat="identity",color=cpalete[1], size=1) +
  geom_point(data= dfplt.month,  aes(x=date, y=RInt, color="Data"), size =2.8) +
  labs(x="", y="\n Reported travel-related cases (monthly)", title="International" ) +
  ylim(0, 30) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11))+
  labs(x="", y="\n Reported travel-related cases (monthly)",
       title="Comparison of model predictions to reported travel-related cases \n International") +
  scale_fill_manual(name = "", labels = c("Best model","Mechanistic model", "Data"), values = c( cpalete[8],cpalete[7], cpalete[1])) +
  scale_color_manual(values = plt1)+
  guides(color = guide_legend(title = ""), size=rel(1.4)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)),
                          legend.position="bottom",
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))


# 3.4. monthly Canada
df2.Dom <- select(dfplt.month, c(date, bCA, CA))
df2.Dom <- melt(df2.Dom, id=c("date"))

p4 <- ggplot() +
  geom_bar(data=df2.Dom, aes(x=date, y=value, fill=variable), stat="identity", position=position_dodge(),alpha=.8 ) +
  geom_line(data= dfplt.month, aes(x=date, y= RCA),stat="identity",color=cpalete[2], size=1) +
  geom_point(data= dfplt.month,  aes(x=date, y= RCA, color="Data"), size =2.8) +
  ylim(0, 120) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11))+
  labs(x="", y="\n Reported travel-related cases (monthly)",
       title="Comparison of model predictions to reported travel-related cases \n Canada") +
  scale_color_manual(values = plt2)+
  scale_fill_manual(name = "", labels = c( "Best model","Mechanistic model", "Data"), values = c( cpalete[8],cpalete[7], cpalete[2])) +
  guides(color = guide_legend(title = ""), size=rel(1.4)) +
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)),
                          legend.position="bottom",
                          legend.text=element_text(size=rel(1.2)),
                          plot.title=element_text(size=rel(1.3), face="bold"),
                          axis.title = element_text(size=rel(1.1)),
                          axis.text.y = element_text(size=rel(1.2)))



p2+p1 + plot_annotation(tag_levels = 'A')
ggsave("Figure/figureApp1.png", width = 15, , height = 6, dpi =500)


(m1+m2)/(p4+p3) + plot_annotation(tag_levels = 'A')
ggsave("Figure/Figure3.png", width = 15, , height = 10, dpi =500)


