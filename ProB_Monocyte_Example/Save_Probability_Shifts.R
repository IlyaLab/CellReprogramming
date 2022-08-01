library(stringr)
library(GA)
intervention_folder="Intervention_State_Probabilities"
nbits=30
nsamples=500

departure_atts=read.table("departure_atts.csv", sep=",", row.names = NULL, header=TRUE)
departure_atts=departure_atts[,2:ncol(departure_atts)]
destination_atts=read.table("destination_atts.csv", sep=",", row.names = NULL, header=TRUE)
destination_atts=destination_atts[,2:ncol(destination_atts)]

original_atts=cbind(apply(departure_atts, 1,paste, collapse=""), apply(destination_atts, 1,paste, collapse=""))
colnames(original_atts)=c("Departure", "Destination")

found_attractor_folder="State_Probabilities"

shift_df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(shift_df) <-c( "Sample", "TF", "Intervention", "Shift")

for (i in 1:nsamples){
    source=original_atts[i,"Departure"]
    target=original_atts[i,"Destination"]
    found_states=read.csv(file=paste(found_attractor_folder, "/",   i ,  "_states.csv", sep=""), colClasses=c("character", "numeric"))                                                                                      
    tot_visits=sum(found_states[,2])
    inds_atts=order(unlist(found_states[,2]), decreasing=TRUE)[1:2]
    rel_mat=found_states[inds_atts, ]
    attractors=unlist(rel_mat[,1])
    source_ini_prob=rel_mat[which(attractors==source),2]/tot_visits
    target_ini_prob=rel_mat[which(attractors==target),2]/tot_visits

  for (bit in 1:nbits){
      new_source <- str_c(str_sub(source, 1, (bit-1)), str_sub(source, bit+1, nchar(source)))
      new_target <- str_c(str_sub(target, 1, (bit-1)), str_sub(target, bit+1, nchar(source)))
      for  (intervention in 0:1){
        interv=read.csv2(paste(intervention_folder,  "/", i, "_", bit, "_", intervention, "_states.csv", sep=""),colClasses=c("character", "numeric"),  sep=",")
        visit_sum=sum(interv[,"second"])
        source_ind=which(interv[,1]==new_source)
        if (length(source_ind)>0){ 
          source_after_prob=interv[source_ind,2]/visit_sum
        }
        else{
          source_after_prob=0
        }
        target_ind=which(interv[,1]==new_target)
        if ((length(target_ind)>0)){ 
          target_after_prob=interv[target_ind,2]/visit_sum
        }
        else{
          target_after_prob=0
        }

        shift=((source_ini_prob-source_after_prob) + (target_after_prob- target_ini_prob))/2
        shift_df[nrow(shift_df) + 1,] <-c(i, bit, intervention, shift)
      }
    }
}

write.csv(shift_df, "prob_shifts_after_intervention.csv",row.names = FALSE)
