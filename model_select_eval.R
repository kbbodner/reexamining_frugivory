# Packages:
library(here) # paths to project's files
library(dplyr) # dataframes 
library(rlist) # lists
library(glmmTMB) # model fitting
library(DHARMa) # model evaluation/diagnosis plot
library(MuMIn) # model selection
library(parallel) # parallel functions
library(pryr) # editing formulas
library(stringr) # editing strings
library(performance)

# Organizing paths to files:
wd <- list() # working directory
wd$data <- here("data") # data directory

# Loading files:

#  Only species with gape size data:
load(file=(paste(sep="/", wd$data, "dataframe_gape_scaled.RData"))) 

# --------------------------------------------------------------------------------------------------------- #
source('00_functionS.R') # calling out functions

# specify the number of cores
ncores <- detectCores()


################################################
#AUTHORS "BEST-PERFORMING" MODEL (ACCORDING TO AIC) FROM DREDGE FUNCTION

# Create the Authors "best" full model with all appropriate terms
best_full_model_gape <- glmmTMB(
  int_frequency ~ (1|net_id/bird_species) + (1|bird_species) + # random effects
    mismatching_gape + dist_edge + degree_frugivory + migratory_status +
    mismatching_gape*dist_edge +  mismatching_gape*degree_frugivory +
    mismatching_gape*migratory_status + dist_edge*degree_frugivory +
    dist_edge*migratory_status + mismatching_gape*dist_edge*degree_frugivory +# fixed effects
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
# Dispersion test:
testDispersion(simulation_best_full_model)
# Plot of the residuals against the predicted value: 
plotResiduals(simulation_best_full_model)
# QQplot:
plotQQunif(simulation_best_full_model,
           testUniformity = T,
           testOutliers = T,
           testDispersion = T)

#testOutliers(simulation_best_full_model, type="bootstrap")


#VARIANCE EXPLAINED

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


#########################################################
# OFFSET-ONLY MODEL

#offset and random effects structure only model
offset_model_gape <- glmmTMB(
  int_frequency ~ (1|net_id/bird_species) + (1|bird_species) + # random effects
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


simulation_best_model_offset <- simulateResiduals(offset_model_gape)
# Dispersion test:
testDispersion(simulation_best_model_offset)
# Plot of the residuals against the predicted value: 
plotResiduals(simulation_best_model_offset)
# QQplot:
plotQQunif(simulation_best_model_offset,
           testUniformity = T,
           testOutliers = T,
           testDispersion = T)

#testOutliers(simulation_best_model_offset, type="bootstrap")



#Variance explained
# Marginal R2 (fixed-effects only):
offset_predict_gape_marginal <- predict(offset_model_gape,
                                      re.form = NA, # population-level predictions (i.e., setting all random effects to zero)
                                      type="response")

cor_marginal <- cor(offset_predict_gape_marginal, offset_model_gape$frame$int_frequency)
cor_marginal^2

# Conditional R2 (fixed and random effects):
offset_predict_gape_conditional <- predict(offset_model_gape,
                                         re.form = NULL, # individual-level predictions
                                         type="response")
cor_conditional <- cor(offset_predict_gape_conditional, offset_model_gape$frame$int_frequency)
cor_conditional^2

# Same result is obtained using the r2_zeroinflated function:
r2_zeroinflated(offset_model_gape, method="correlation")

#############################################################

#AIC COMPARISON
AIC(best_full_model_gape,offset_model_gape)

# Calculate BIC for both models
bic_full_model <- BIC(best_full_model_gape)
bic_offset_model <- BIC(offset_model_gape)

# Print BIC values
cat("BIC for Full Model:", bic_full_model, "\n")
cat("BIC for Offset Model:", bic_offset_model, "\n")

# Calculate Delta BIC
delta_bic_full_vs_offset <- bic_full_model - bic_offset_model
cat("Delta BIC (Full Model vs Offset Model):", delta_bic_full_vs_offset, "\n")











