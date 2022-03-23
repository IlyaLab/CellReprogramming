
library("TripleR")
library("RobustRankAggreg")
n_sets=500
input_folder="Output"

setwd(input_folder)
orig_atts=read.table("Attractors.csv", sep=",", header=FALSE)
attractors_found=matrix(0,n_sets,2)
attractor_probabilities=matrix(0,n_sets,2)
for (i in 1:n_sets){ 
  attractors_found[i,]=as.numeric(read.table(paste(i-1,"/", "found_attractors.csv", sep=""), sep=",", header=FALSE))
  attractor_probabilities[i,]=as.numeric(read.table(paste(i-1, "/",  "attractor_probability.csv", sep=""), sep=",", header=FALSE))         
}
hist(as.vector(attractor_probabilities), main = "" , xlab="Cell Type Probabilities")
#######-------#####-------

TFs=unlist(read.table("TF_list.csv"))

results=list()
melted_results=list()

for (i in 1:n_sets){ 
  results[[i]]=read.table(paste(i-1, "/intervention_result.csv",sep=""), sep=",",  header=FALSE)
  colnames(results[[i]])=c("OFF", "ON")
  rownames(results[[i]])=TFs
  temp=matrix2long(results[[i]], new.ids=FALSE, var.id="SHIFT")
  colnames(temp)=c("TF" , "INTERV", "SHIFT")
  temp2=cbind(paste(temp$TF,temp$INTERV), temp$SHIFT)
  colnames(temp2)=c("TF_INTV", "SHIFT")
  temp2=temp2[order(-as.numeric(as.vector(temp2[,"SHIFT"]))), ]
  melted_results[[i]]=temp2[,1]
}

aggreg_list=aggregateRanks(melted_results, full = TRUE)

-1*log(aggreg_list[1:4,2])
aggreg_list[1:4,1]

counts <- -1*log(aggreg_list[1:4,2])
barplot(counts, main="Intervention Effect on Transdifferentiation", ylab="Weight",xlab="Intervention" , names.arg=aggreg_list[1:4,1], cex.names=0.6)

ebf1_down=cebpb_up=tcf3_down=stat3_up=vector(length=n_sets)

for (i in 1:n_sets){
    ebf1_down[i]=results[[i]]["EBF1", "OFF"]
    cebpb_up[i]=results[[i]]["CEBPB", "ON"]
    tcf3_down[i]=results[[i]]["TCF3", "OFF"]
    stat3_up[i]=results[[i]]["STAT3", "ON"]
}

ebf1_down_df=data.frame(rep.int("EBF1 OFF", n_sets), ebf1_down)
cebpb_up_df=data.frame(rep.int("CEBPB ON", n_sets), cebpb_up)
tcf3_down_df=data.frame(rep.int("TCF3 OFF", n_sets), tcf3_down)
stat3_up_df=data.frame(rep.int("STAT3 ON", n_sets), stat3_up)
names(ebf1_down_df)=names(cebpb_up_df)=names(tcf3_down_df)=names(stat3_up_df)=c("Intervention", "Shift")
all_res=data.frame(rbind(ebf1_down_df, cebpb_up_df, tcf3_down_df, stat3_up_df))

x1 <- all_res$Shift[all_res$Intervention=="EBF1 OFF"]
x2 <- all_res$Shift[all_res$Intervention=="CEBPB ON"]
x3 <- all_res$Shift[all_res$Intervention=="TCF3 OFF"]
x4 <- all_res$Shift[all_res$Intervention=="STAT3 ON"]

boxp=ggplot(all_res, aes(x=reorder(Intervention,-Shift), y=Shift) )+ 
labs(y="Probability Shift (PS)", x="Intervention")+
geom_boxplot() + theme_bw()

ggsave(boxp, file="shift_boxplot.png", width = 7, height = 4, units ="in",  dpi = 300)
ggsave(boxp, file="shift_boxplot.pdf", width = 7, height = 4, units ="in",  dpi = 300)

