if (input_matrix != "") {
  print("Reading filtered matrices per subset")
  load(input_matrix)
  if (module_comparisons=="YES") {
    if (method=="Shallow"){
      number_of_subsets=length(subset_ids)
      for (i in 1:number_of_subsets){
        filtered_trait=filtered_trait_objs[[i]]
        assign(paste('filtered_trait',subset_ids[i],sep='_'), filtered_trait)
        print(subset_ids[i])
      } 
    } else {
      j=1
      for (i in (1+number_of_subsets):length(filtered_trait_objs)){
        filtered_trait=filtered_trait_objs[[i]]
        assign(paste('filtered_trait',subset_ids[j],sep='_'), filtered_trait)
        j=j+1
        print(subset_ids[i])
      }
    }
  }
} else {
  stop("No matrices provided")
}