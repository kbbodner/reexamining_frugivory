# Packages:
library(here) # paths to project's files
library(glmmTMB) # model fitting
library(DHARMa) # model evaluation/diagnosis plot
library(MuMIn) # model selection
library(parallel) # parallel functions
library(performance)

# Organizing paths to files:
wd <- list() # working directory
wd$data <- here("data") # data directory

# Loading files:

#  Only species with gape size data:
load(file=(paste(sep="/", wd$data, "dataframe_gape_scaled.RData"))) 

# specify the number of cores
ncores <- detectCores()


################################################
#AUTHORS "BEST-PERFORMING" MODEL (ACCORDING TO AIC) FROM DREDGE FUNCTION

#extra fixed effect that is in equally good model: dist_edge*migratory_status according to AIC

# Create the Authors "best" full model with all appropriate terms
best_full_model_gape <- glmmTMB(
  int_frequency ~ (1|net_id/bird_species) + (1|bird_species) + # random effects
    mismatching_gape + dist_edge + degree_frugivory + migratory_status +
    mismatching_gape*dist_edge +  mismatching_gape*degree_frugivory +
    mismatching_gape*migratory_status + dist_edge*degree_frugivory +
     mismatching_gape*dist_edge*degree_frugivory +# fixed effects
    offset(log(local_plant_availability)), # offset (plant availability)
  dispformula = ~ 1, # dispersion formula
  ziformula = ~ 1 + (1|net_id), # zi formula
  family = nbinom2, # negative binomial (used to be nbinom2)
  data=dataframe_gape_scaled, # dataset
  verbose=F,
  se=TRUE,
  control = glmmTMBControl(
    optCtrl = list(iter.max=1000, eval.max=1000),
    profile=TRUE,
    parallel=ncores
  )
)

null_model_gape <- glmmTMB(
  int_frequency ~ (1|net_id/bird_species) +  (1|bird_species) +# random effects
  offset(log(local_plant_availability)), # offset (plant availability)
  dispformula = ~ 1, # dispersion formula
  ziformula = ~ 1 + (1|net_id), # zi formula
  family = nbinom2, # negative binomial (used to be nbinom2)
  data=dataframe_gape_scaled, # dataset
  verbose=F,
  se=TRUE,
  control = glmmTMBControl(
    optCtrl = list(iter.max=1000, eval.max=1000),
    profile=TRUE,
    parallel=ncores
  )
)

#TESTING MODEL ASSUMPTIONS

simulation_best_full_model <- simulateResiduals(best_full_model_gape)

####### Individual tests #######

# KS test for correct distribution of residuals
testUniformity(simulation_best_full_model)

#Testing for Outliers
testOutliers(simulation_best_full_model, type="bootstrap")

# KS test for correct distribution within and between groups
testCategorical(simulation_best_full_model, dataframe_gape_scaled$net_id)

# Dispersion test - for details see ?testDispersion
testDispersion(simulation_best_full_model) # tests under and overdispersion

# testing zero inflation
testZeroInflation(simulation_best_full_model)

# Plot of the residuals against the predicted value: 
plotResiduals(simulation_best_full_model)

# QQplot:
plotQQunif(simulation_best_full_model,
           testUniformity = T,
           testOutliers = T,
           testDispersion = T)


######## VARIANCE EXPLAINED (BY MARTINS ET AL) ##########

# Marginal R2 (fixed-effects only):
best_full_predict_gape_marginal <- predict(best_full_model_gape,
                                      re.form = NA, # population-level predictions (i.e., setting all random effects to zero)
                                      type="response")

best_full_cor_marginal <- cor(best_full_predict_gape_marginal, best_full_model_gape$frame$int_frequency)
best_full_cor_marginal^2

# Conditional R2 (fixed and random effects):
best_full_predict_gape_conditional <- predict(best_full_model_gape,
                                         re.form = NULL, # individual-level predictions
                                         type="response")
best_full_cor_conditional <- cor(best_full_predict_gape_conditional, best_full_model_gape$frame$int_frequency)
best_full_cor_conditional^2

# Same result is obtained using the r2_zeroinflated function:
r2_zeroinflated(best_full_model_gape, method="correlation")

#######

# Marginal R2 (fixed-effects only):
null_model_gape_marginal <- predict(null_model_gape,
                                           re.form = NA, # population-level predictions (i.e., setting all random effects to zero)
                                           type="response")

null_cor_marginal <- cor(null_model_gape_marginal, null_model_gape$frame$int_frequency)
null_cor_marginal^2

# Conditional R2 (fixed and random effects):
null_predict_gape_conditional <- predict(null_model_gape,
                                              re.form = NULL, # individual-level predictions
                                              type="response")
null_cor_conditional <- cor(null_predict_gape_conditional, null_model_gape$frame$int_frequency)
null_cor_conditional^2

# Same result is obtained using the r2_zeroinflated function:
r2_zeroinflated(null_model_gape, method="correlation")

#############################################################

#AIC COMPARISON
AIC(best_full_model_gape, null_model_gape)
AIC(best_full_model_gape) - AIC(null_model_gape)

# Calculate BIC for both models
BIC(best_full_model_gape,null_model_gape)

bic_best_full_model_gape <- BIC(best_full_model_gape)
bic_null_model_gape <- BIC(null_model_gape)

# Calculate Delta BIC
delta_bic_full_vs_nested <- bic_best_full_model_gape - bic_null_model_gape

#print
delta_bic_full_vs_nested

#############################################################

#MEAN ABSOLUTE ERROR

# Make predictions on the dataset
predictions_best_full_model_gape <- predict(best_full_model_gape, type = "response")
predictions_null_model_gape <- predict(null_model_gape, type = "response")

#MEAN ABSOLUTE ERROR

# Calculate the absolute errors
absolute_errors_best_full_model_gape <- abs(dataframe_gape_scaled$int_frequency - predictions_best_full_model_gape)
absolute_errors_null_model_gape <- abs(dataframe_gape_scaled$int_frequency - predictions_null_model_gape)

# Compute the Mean Absolute Error
mae_best_full_model_gape <- mean(absolute_errors_best_full_model_gape)
mae_null_model_gape <- mean(absolute_errors_null_model_gape)

#print
mae_best_full_model_gape
mae_null_model_gape










