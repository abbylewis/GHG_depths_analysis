#' Run one permutation of random forest
#'
#' @param rep_n Which rep are we on
#' @param for_rf_surf Data to use
#' @param num_vars Which variables are numeric (for partial dependence)
#' @param cat_vars Which variables are categorical (for partial dependence)
#'
#' @returns list with model R2, data frame of partial regression results, and data frame of variable importance
#' 
permutation_rf <- function(rep_n, for_rf_surf, num_vars, cat_vars){
  for_rf_surf_lim <- for_rf_surf %>%
    group_by(LakeID) %>%
    slice_sample(n = 1) %>%
    ungroup() %>%
    select(-LakeID)
  rf_surf <- randomForest(value ~ ., data = for_rf_surf_lim, 
                          proximity=TRUE, importance = T) 
  surf_r2 = rf_surf$rsq[length(rf_surf$rsq)]
  
  # Calculate partial correlations
  partials_df_surf <- pdp::partial(rf_surf, pred.var = num_vars[1], 
                                   train = as.data.frame(for_rf_surf)) %>%
    pivot_longer(-yhat, names_to = "var", values_to = "x") %>%
    rename(y = yhat)
  for(i in 2:length(num_vars)){
    new <- pdp::partial(rf_surf, pred.var = num_vars[i], 
                        train = as.data.frame(for_rf_surf)) %>%
      pivot_longer(-yhat, names_to = "var", values_to = "x") %>%
      rename(y = yhat)
    partials_df_surf = rbind(partials_df_surf, new)
  }
  partials_df_surf <- partials_df_surf %>%
    mutate(rep = rep_n)
  
  ImpData_surf <- as.data.frame(randomForest::importance(rf_surf)) %>% 
    mutate(rep = rep_n)
  ImpData_surf$Var.Names <- row.names(ImpData_surf)
  
  return(list(surf_r2, partials_df_surf, ImpData_surf))
}