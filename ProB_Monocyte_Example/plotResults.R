library("TripleR")
library("RobustRankAggreg")
library("dplyr")
n_sets=500

TFs=unlist(read.table("TF_list.csv"))
names(TFs)=NULL
results=list()
melted_results=list()
shift_df=read.csv("prob_shifts_after_intervention.csv")

#"Sample"  "TF" "Intervention" "Shift"  

shift_df$TF=TFs[shift_df$TF]

shift_df <- shift_df %>% 
  mutate("Intervention" = if_else( Intervention==0, "OFF", "ON"))


for (i in 1:n_sets){ 
    temp=subset(shift_df, Sample==i, select=c("TF", "Intervention", "Shift"))
    colnames(temp)=c("TF" , "INTERV", "SHIFT")
    temp2=cbind(paste(temp$TF,temp$INTERV), temp$SHIFT)
    colnames(temp2)=c("TF_INTV", "SHIFT")
    temp2=temp2[order(-as.numeric(as.vector(temp2[,"SHIFT"]))), ]
    melted_results[[i]]=temp2[,1]
  }
aggreg_list=aggregateRanks(melted_results, full = TRUE)

Interventions=list()
x=sapply(aggreg_list[1:4, "Name"], function(x){str_split(x," ")})
hm=0
for (s in 1:4){
  hm=hm+1
  shift=subset(shift_df, (TF==x[[s]][1] &  Intervention==x[[s]][2]), select=c("Shift"))
  xs=data.frame( "Intervention"=aggreg_list[s, "Name"], "Shift"=shift)
  Interventions[[s]] <- xs   
}
plot_data=do.call(rbind, Interventions)

boxp=ggplot(plot_data, aes(x=reorder(Intervention,-Shift), y=Shift) )+ 
  labs(y="Probability Shift (PS)", x="Intervention")+
  geom_boxplot() + theme_bw()
ggsave(boxp, file="shift_boxplot.png", width = 7, height = 4, units ="in",  dpi = 300)

