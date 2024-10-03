
###########################################################################

#                                 FUNCTIONS                               #

###########################################################################



#### ------------------ Network data --------------- ####

# Function to load plant-frugivore networks:
upload_networks <- function (directory_name)   {
  
  load(file=(paste(sep="/", directory_name, "mynetworks.RData")))
  load(file=(paste(sep="/", wd$networks, "list_networks.RData")))
  
  original_networks <- mynetworks[list_networks]
  original_networks <- original_networks[order(names(original_networks))]
  
  return(original_networks)
}

# Function to remove 'problematic' species from the networks:
removal_problematic <- function (networks){
  load(file=(paste(sep="/", wd$networks, "problematic_species.RData")))
  
  updated_networks<- networks
  
  for(i in 1:length(networks)){
    remove_positions <-  match(problematic_species, row.names(networks[[i]]))
    remove_positions <- remove_positions[!is.na(remove_positions)]
    if(length(remove_positions)==0){
      updated_networks[[i]] <- networks[[i]]
    }else{
      removal_plant <- networks[[i]][-remove_positions,]
      removal_bird <- which(colSums(removal_plant)==0)
      if(length(removal_bird)==0){
        updated_networks[[i]] <- networks[[i]][-remove_positions,]
      } else{
        updated_networks[[i]] <- networks[[i]][-remove_positions, -removal_bird]
      }
    }
  }
  return(updated_networks)
}

# Function to update bird species names:
update_birds <- function(networks){
  
  load( file=(paste(sep="/", wd$networks, "bird_tax_former.RData")))
  load(file=(paste(sep="/", wd$networks, "bird_tax_new.RData")))
  
  updated_networks <- networks
  net_upd_birds <- names(bird_tax_new)
  for(k in net_upd_birds){
    positions_bird <- match(bird_tax_former[[k]], colnames(networks[[k]]))
    
    if(length(positions_bird)==0){
      colnames(updated_networks[[k]]) <- colnames(networks[[k]])
    } else{
      colnames(updated_networks[[k]])[positions_bird] <- bird_tax_new[[k]]
    }
  }
  return(updated_networks)
}

# Function to update plant species names:
update_plants <- function(networks){
  
  load(file=(paste(sep="/", wd$networks, "plant_tax_former.RData")))
  load(file=(paste(sep="/", wd$networks, "plant_tax_new.RData")))
  
  updated_networks <- networks
  net_upd_plants <- names(plant_tax_new)
  for(k in net_upd_plants){
    positions_plant <- match(plant_tax_former[[k]], row.names(networks[[k]]))
    
    if(length(positions_plant)==0){
      row.names(updated_networks[[k]]) <- row.names(networks[[k]])
    } else{
      row.names(updated_networks[[k]])[positions_plant] <- plant_tax_new[[k]]
    }
  }
  return(updated_networks)
}


#### ------------------ Model selection ----------- ####

# Function to run model selection:
model_selection <- function(main_combinations, combinations_moderators, dataframe){
  model_selection_main <- list()
  model_selection_mod <- list()
  
  #### Combinations - focal effects:
  for(i in 1:length(main_combinations)){
    
    # Fitting model:
    sub_model <- eval(main_combinations[[i]]) 
    
    coeffts1 <- coeffs(sub_model)
    AIC_model1 <- AIC(sub_model)
    AICc_model1 <- MuMIn::AICc(sub_model)
    loglik1 <- logLik(sub_model)  
    
    # Removing object to save space/memory
    rm(sub_model)
    
    model_selection_main [[i]] <- data.frame(coefs = coeffts1, AIC = AIC_model1,AICc =AICc_model1, LogLik = loglik1)  
    
    # Progress bar:
    Sys.sleep(0.01)
    
    # Print progress
    print(paste("Finished", i, "of", length(main_combinations), "(main variables)"))
    
  }
  
  #### Combinations - moderating effects:
  re <- "(1|net_id/bird_species) + (1|bird_species)" # random effects
  offset_term <- "offset(log(local_plant_availability))" # offset
  
  for(k in 1:length(combinations_moderators)){
    formula_mod <- as.formula(paste0('int_frequency~', combinations_moderators[[k]], 
                                     '+', re,'+', offset_term ))
    
    # Fitting model:
    model_moderators <- glmmTMB(formula_mod, # formula
                                dispformula = ~ 1, # dispersion formula
                                ziformula = ~ 1 + (1|net_id), # zi formula
                                family = nbinom2, # negative binomial
                                data=dataframe, # dataset
                                verbose=FALSE,
                                se=TRUE,
                                control = glmmTMBControl(
                                  optCtrl = list(iter.max=1000, eval.max=1000),
                                  profile=TRUE,
                                  parallel=ncores
                                ))
    
    coeffts2 <- coeffs(model_moderators)
    AIC_model2 <- AIC(model_moderators)
    AICc_model2 <- MuMIn::AICc(model_moderators)
    loglik2 <- logLik(model_moderators)  
    
    rm(model_moderators)
    
    model_selection_mod [[k]] <- data.frame(coefs = coeffts2, AIC = AIC_model2,AICc =AICc_model2, LogLik = loglik2)  
    
    # Progress bar:
    Sys.sleep(0.01)
    
    # Print progress
    print(paste("Finished", k, "of", length(combinations_moderators), "(with moderators)"))
    
  }
  
  model_selection_list <- c(model_selection_main, model_selection_mod)
  
  names(model_selection_list) <- seq(1,length(model_selection_list))
  
  names_all1 <- unique(unlist(lapply(model_selection_main, row.names)))
  names_all2 <- unique(unlist(lapply(model_selection_mod, row.names)))
  
  names_all <- unique(c(names_all1, names_all2))
  out_names <- unique(unlist(lapply(model_selection_list, colnames)))[-1]
  col_table <- c(names_all, out_names)
  model_selection_table <- matrix(NA, ncol=length(col_table), nrow=length(model_selection_list))
  colnames(model_selection_table) <- col_table
  row.names(model_selection_table) <- names(model_selection_list)
  
  model_selection_table <- as.data.frame(model_selection_table)
  
  for(j in 1:nrow(model_selection_table)){
    positions <- match(row.names(model_selection_list[[j]]), colnames(model_selection_table))
    model_selection_table[j,positions] <- model_selection_list[[j]]$coefs 
    model_selection_table$AIC[j] <- unique(model_selection_list[[j]]$AIC)
    model_selection_table$AICc[j] <- unique(model_selection_list[[j]]$AICc)
    model_selection_table$LogLik[j] <- unique(model_selection_list[[j]]$LogLik)
    
  }
  
  model_selection_table <- model_selection_table%>%arrange(AIC)%>%
    relocate("zi((Int))", .before = AIC)
  
  model_selection_table$deltaAIC <- NA
  
  for(l in 1:nrow(model_selection_table)){
    model_selection_table$deltaAIC[l]<- model_selection_table$AIC[l] -model_selection_table$AIC[1]
  }
  return(model_selection_table)
}


#  Function to organize model selection table:
organize_table <- function(model_table, std=FALSE){
  
  var_subset <- colnames(model_table)
  remove <- c(grep('((Int))', var_subset),
              grep('LogLik', var_subset),
              grep('AIC', var_subset))
  retained <- var_subset[-remove]
  
  retained_variables <- gsub("[\\(\\)]", "", regmatches(retained, gregexpr("\\(.*?\\)", retained)))
  split_var <- function(x){
    str_split(x, ':')
  }
  position_migratory <- grep('migratory_status', retained_variables)
  if(length(position_migratory)>0){
    retained_variables <- str_replace( retained_variables,'migratory_statusnon-migrant', 'migratory_status')
  }
  count_var <- lapply(retained_variables, split_var)
  count_var2 <-lapply(count_var, unlist)
  count_var_final <- lapply(count_var2, sort)
  
  names(count_var_final) <- retained
  
  for(i in 1:length(count_var_final)){
    count_var_final[[i]] <- paste(count_var_final[[i]], collapse=":")
    position <- match(names(count_var_final)[i], colnames(model_table))
    colnames(model_table)[position]<-count_var_final[[i]] 
  }
  
  unique_pred <- unique(colnames(model_table))
  remove2 <- c(grep('((Int))', unique_pred),
               grep('LogLik', unique_pred),
               grep('AIC', unique_pred))
  
  unique_pred <- unique_pred[-remove2]
  
  table_zeroes <-model_table %>% replace(is.na(.), 0)
  colunm_list <- list()
  for(k in 1:length(unique_pred)){
    position <- which(!is.na(match(colnames(model_table),unique_pred[k])))
    if(length(position)>1){
      combined_colunms <- rowSums(table_zeroes[,position])
    } else{
      combined_colunms <- table_zeroes[,position] 
    }
    combined_colunms[combined_colunms == 0] <- NA
    data_table_comb  <- as.data.frame(combined_colunms) 
    colnames(data_table_comb) <- unique_pred[k]
    colunm_list[[k]] <- data_table_comb
  }
  org_table <- list.cbind(colunm_list)
  org_table <- cbind(org_table, model_table[remove])
  if(std==TRUE){
    organized_table <- org_table%>%
      relocate("cond((Int))", .before = "std_dist_edge")
  } else{
    organized_table <- org_table%>%
      relocate("cond((Int))", .before = "dist_edge")
}
  return(organized_table)
  
}


# Function to retrieve the retained predictors from best-fitting model(s):
retained_predictors <- function(model_selection_table){
  best_fitting_subset <- filter(model_selection_table, deltaAIC<2)
  best_fitting_subset[is.na(best_fitting_subset)] <- 0
  var_subset <- colnames(best_fitting_subset)[!(colSums(best_fitting_subset)==0)]
  remove <- c(grep('((Int))', var_subset),
              grep('LogLik', var_subset),
              grep('AIC', var_subset))
  retained <- var_subset[-remove]
  
  retained_variables <- gsub("[\\(\\)]", "", regmatches(retained, gregexpr("\\(.*?\\)", retained)))
  split_var <- function(x){
    str_split(x, ':')
  }
  position_migratory <- grep('migratory_status', retained_variables)
  if(length(position_migratory)>0){
    retained_variables <- str_replace( retained_variables,'migratory_statusnon-migrant', 'migratory_status')
  }
  
  count_var <- lapply(retained_variables, split_var)
  position_3 <- which(lapply(lapply(count_var, '[[', 1), length)>2)
  
  for(i in position_3){
    retained_variables[i] <- paste0('(', str_replace_all(retained_variables[i], ':', '+'), 
                                    ')^3')
  }
  return(retained_variables)
}

# Function to retrieve the top best-fitting model (i.e., lowest AIC):
get_best_model <- function (combined_model,model_dredge_table){
  model_dredge_comb <- dredge(global.model = combined_model, 
                              fixed = "cond(offset(log(local_plant_availability)))",
                              evaluate = F)
  
  best_fit_formula <- model_dredge_comb[grep(row.names(model_dredge_table)[1], names(model_dredge_comb))]
  best_fit_model <-  eval(best_fit_formula[[1]]) 
  return(best_fit_model)
}

# Function to collapse formula:
formula_collapse <- function(x){
  paste0(x, collapse='+')
}

#### ------------------ Sensitivity analyses ----------- ####

# Function to generate dataframes with different penalizations for trait mismatching:
get_dataframes_sensitivity <- function(dataframe, penalty_number, one.sided=T){
  dataframe_list <- list()
  
  dataframe_pos_matching <- dataframe
  dataframe_pos_matching$mismatching_gape_raw <- dataframe_pos_matching$gape_size-dataframe_pos_matching$fruit_diameter
  dataframe_pos_matching$mismatching_gape_original <- dataframe_pos_matching$mismatching_gape
  dataframe_list[[1]] <- dataframe_pos_matching %>% 
    mutate_at(c('mismatching_gape','dist_edge', 'dist_elevational_limit', 
                'degree_frugivory', 'human_footprint', 'abs_lat' ), scale)%>%
    select(net_id, net_lat, net_lon, locality, bird_species, gape_size,
           plant_species, fruit_diameter, local_plant_availability, dist_edge, dist_elevational_limit,
           degree_frugivory, migratory_status, human_footprint, abs_lat, mismatching_gape_raw,
           mismatching_gape_original, mismatching_gape, int_frequency)
  
  for(i in 1:penalty_number){
    new_dataframe <- dataframe_pos_matching
    
    new_dataframe$mismatching_gape_original <- dataframe_pos_matching$mismatching_gape
    new_dataframe$penalty <- i
    
    if(one.sided==T){ # one-sided analysis (i.e., change in peak)
      new_dataframe$mismatching_gape <- abs(dataframe_pos_matching$mismatching_gape_raw-i)
      new_dataframe$mismatching_gape_non_scaled <- abs(dataframe_pos_matching$mismatching_gape_raw-i)
      
    }else{ # two-sided analysis (i.e., increase in range)
      new_zeroes <- which(dataframe_pos_matching$mismatching_gape<=i)
      new_dataframe$mismatching_gape[new_zeroes]<- 0
      
      position_out_range <- which(dataframe_pos_matching$mismatching_gape>i)
      new_dataframe$mismatching_gape[position_out_range]<- dataframe_pos_matching$mismatching_gape[position_out_range]-i 
      new_dataframe$mismatching_gape_non_scaled <- abs(new_dataframe$mismatching_gape)
    }
    dataframe_list[[i+1]] <- new_dataframe %>% 
      mutate_at(c('mismatching_gape','dist_edge', 'dist_elevational_limit', 
                  'degree_frugivory', 'human_footprint', 'abs_lat' ), scale)%>%
      select(net_id, net_lat, net_lon, locality, bird_species, gape_size,
             plant_species, fruit_diameter, local_plant_availability, dist_edge, dist_elevational_limit,
             degree_frugivory, migratory_status, human_footprint, abs_lat, mismatching_gape_raw, mismatching_gape_original,
             penalty,mismatching_gape_original, mismatching_gape_non_scaled, mismatching_gape,int_frequency)
  }
  return(dataframe_list)
  
}

# Function for running sensitivity analyses:
run_sensitivity_analysis <- function(best_model, dataframe_list, cores){
  model_sensitivity <- list()
  summary_sensitivity <- list()
  coeffts <- list()
  output_table <- list()
  
  formula_sensitivity <- best_model$call
  
  for(i in 1:length(dataframe_list)){
    dataframe_gape_scaled <- dataframe_list[[i]]
    
    model_sensitivity[[i]]  <- eval(formula_sensitivity)
    
    # Saving in the 'sensitivity' folder (optional)
    #save(model_sensitivity, file=paste(sep='/', wd$sensitivity, 'model_sensitivity.RData'))
    
    # Model summary:
    summary_sensitivity[[i]] <- summary(model_sensitivity[[i]])
    predictors <- summary_sensitivity[[i]]$coefficients$cond
    coeffts[[i]] <- predictors[(2:nrow(predictors)),1:2]
    
    output_table[[i]] <-   as.data.frame(cbind(variable = row.names(coeffts[[i]]), penalty = i-1, coeffts[[i]]))
    
    # Saving in the 'sensitivity' folder (optional)
    #save(output_table, file=paste(sep='/', wd$sensitivity, 'output_table.RData'))
    
    # Progress bar:
    Sys.sleep(0.01)
    
    # Print progress
    print(paste("Finished", i, "of", length(dataframe_list)))
    
  }
  outuput_sensitivity <- data.table::rbindlist(output_table)
  
  return(outuput_sensitivity)
}

# Function to get the statistics for each penalty dataframe:
get_statistics <- function(table, penalty_number){
  coef_names <- unique(table$variable)
  
  # Calculating confidence intervals:
  table$lower_CI <- NA  
  table$upper_CI <-NA  
  
  for(i in 1:nrow(table)){
    estimate <- as.numeric(table$Estimate[i])
    std_error <- as.numeric(table$`Std. Error`[i])
    table$lower_CI[i] <- estimate-1.96*std_error  
    table$upper_CI[i] <- estimate+1.96*std_error  
  }
  
  dataframes_penalties <- list()
  
  for(i in 1:(penalty_number+1)){
    dataframes_penalties[[i]] <-   filter(table, penalty==i-1)
  }
  
  stat_list <- list()
  for(k in 1:length(coef_names)){
    stat_list[[k]]   <- data.frame(estimate = as.numeric(unlist(lapply(dataframes_penalties, function(x) x[k,3]), use.names = F)),
                                   std_error = as.numeric(unlist(lapply(dataframes_penalties, function(x) x[k,4]), use.names = F)),
                                   lower_CI = as.numeric(unlist(lapply(dataframes_penalties, function(x) x[k,5]), use.names = F)),
                                   upper_CI = as.numeric(unlist(lapply(dataframes_penalties, function(x) x[k,6]), use.names = F)),
                                   penalty= c(0:(length(dataframes_penalties)-1)),
                                   colour="black")
  }
  names(stat_list) <- c('degree_fugivory',
                        "distance_to_edge",
                        "migratory status",
                        "trait_mismatching",
                        "distance_to_edge:degree_frugivory",
                        "trait_mismatching:degree_frugivory",
                        "trait_mismatching:distance_to_edge",
                        "trait_mismatching:migratory_status",
                        "trait_mismatching:distance_to_edge:degree_frugivory")
  
  return(stat_list)
}


#### ------------------ Plotting ---------------- ####

# Function for generating multiple groups for plotting:
plot_multiple_groups <- function(dataframe, predictor1, predictor2){
  
  number_group <- round(nrow(dataframe)/3)
  group1 <- 1: number_group
  group2 <- (number_group+1):(number_group*2)
  group3 <- ((number_group*2)+1):nrow(dataframe)
  
  # Organizing the dataframe based on the predictor:
  dataframe_arranged1 <- dataframe%>%arrange({{predictor1}})
  
  subset1 <- dataframe_arranged1[group1,] # low
  subset2 <- dataframe_arranged1[group2,] # intermediate
  subset3 <- dataframe_arranged1[group3,] # high
  
  # Organizing based on second predictor:
  subset1 <- arrange(subset1, {{predictor2}}) # low
  subset2 <- arrange(subset2, {{predictor2}}) # intermediate
  subset3 <- arrange(subset3, {{predictor2}}) # high
  
  # Splitting the subgroups into 3-equal sized groups:
  number_subgroup <- round((nrow(subset1)/3))
  group1.1 <- 1: number_subgroup
  group1.2 <- (number_subgroup+1):(number_subgroup*2)
  group1.3 <- ((number_subgroup*2)+1):nrow(subset1)
  
  # Low :
  subset1$subset <- NA
  subset1$subset[group1.1] <- rep('group1', number_subgroup)
  subset1$subset[group1.2] <- rep('group2', number_subgroup)
  subset1$subset[group1.3] <- rep('group3', number_subgroup-1)
  
  # Intermediate
  subset2$subset <- NA
  subset2$subset[group1.1] <- rep('group1', number_subgroup)
  subset2$subset[group1.2] <- rep('group2', number_subgroup)
  subset2$subset[group1.3] <- rep('group3', number_subgroup-1)
  
  # High:
  subset3$subset <- NA
  subset3$subset[group1.1] <- rep('group1', number_subgroup)
  subset3$subset[group1.2] <- rep('group2', number_subgroup)
  subset3$subset[group1.3] <- rep('group3', number_subgroup-1)
  subsets <- list()
  subsets[[1]] <- subset1
  subsets[[2]] <- subset2
  subsets[[3]] <- subset3
  names(subsets) <- c('Low', 'Intermediate', 'High')
  
  return(subsets)
  
}  

# Function for generating groups for plotting:
plot_group <- function(dataframe, predictor){
  
  # Splitting the dataset into 3-equal sized groups:
  number_group <- round((nrow(dataframe)/3))
  group1 <- 1: number_group
  group2 <- (number_group+1):(number_group*2)
  group3 <- ((number_group*2)+1):nrow(dataframe)
  
  # Organizing the dataframe based on the predictor:
  dataframe_arranged1 <- dataframe%>%arrange({{predictor}})
  
  subset1 <- dataframe_arranged1[group1,] # low
  subset2 <- dataframe_arranged1[group2,] # intermediate
  subset3 <- dataframe_arranged1[group3,] # high
  
  subset <- rbind(subset1,subset2,subset3)
  subset$subset <- NA
  subset$subset[group1] <- rep('group1', number_group)
  subset$subset[group2] <- rep('group2', number_group)
  subset$subset[group3] <- rep('group3', length(group3))
  
  return(subset)
  
}  

# Function for adding inset plots into each facet:
inset_density_plot <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) 
{
  layer(data = data, stat = StatIdentity, position = PositionIdentity, 
        geom = ggplot2:::GeomCustomAnn,
        inherit.aes = FALSE, params = list(grob = grob, 
                                           xmin = xmin, xmax = xmax, 
                                           ymin = ymin, ymax = ymax))
}

# Function to get the data necessary for figures S16 and S18
get_sensitivity_data <- function(dataframe_list, penalty_values){
  position_list <- penalty_values+1
  output_list <- list()
  for(i in position_list){
    penalty_dataframe <- data.frame(penalty = dataframe_list[[i]]$penalty,
                                    trait_matching_original = dataframe_list[[i]]$mismatching_gape_original,
                                    trait_matching_penalized = dataframe_list[[i]]$mismatching_gape_non_scaled,
                                    interaction_frequency =  dataframe_list[[i]]$int_frequency)
    zero_frequency <- which(penalty_dataframe$interaction_frequency==0)
    
    penalty_dataframe$colour_var <- 'interaction'
    penalty_dataframe$colour_var[zero_frequency] <- 'no interaction'
    
    output_list[[i]] <- penalty_dataframe
  }
  output <- output_list[position_list]
  names(output) <- paste0("penalty", penalty_values)
  return(output)
}