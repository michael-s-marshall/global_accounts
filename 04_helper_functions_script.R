# helper functions -----------------------------------------

give_vpcs <- function(mlm, group){
  
  as_tibble(summary(mlm)$varcor) %>%
    select(grp,vcov) %>% 
    pivot_wider(names_from = "grp", 
                values_from = "vcov") %>% 
    mutate("{{group}}_vpc" := {{group}} / ({{group}} + Residual),
           unit_vpc = 1 - ({{group}} / ({{group}} + Residual)))
  
}

## [rescale01] Function to rescale a variable from 0 to 1
rescale01 <- function(x, ...) {
  (x - min(x, ...)) / ((max(x, ...)) - min(x, ...))
}

mlm_diagnostics <- function(mlm_mod){
  
  require(patchwork)
  
  diag_plots <- sjPlot::plot_model(mlm_mod, type = "diag")
  
  (diag_plots[[1]] / diag_plots[[3]]) | diag_plots[[4]]
  
}

# my summarise

my_summarise <- function(df, vars){
  
  df %>% 
    summarise(
      across({{vars}}, mean, na.rm = T, .names = "mean_{.col}"),
      across({{vars}}, median, na.rm = T, .names = "median_{.col}"),
      across({{vars}}, sd, na.rm = T, .names = "sd_{.col}")
    )
  
  
}

# group mean centering

group_centre <- function(df, group_var, centre_vars){
  
  group_mc <- function(x){
    x - mean(x, na.rm = T)
  }
  
  df %>% 
    group_by({{group_var}}) %>% 
    mutate(
      across({{centre_vars}},
             group_mc, .names = "{.col}_gc")
    ) %>% 
    ungroup()
  
}

group_mean <- function(df, group_var, mean_vars){
  
  df %>% 
    group_by({{group_var}}) %>% 
    mutate(
      across({{mean_vars}},
             mean, na.rm = T, .names = "{.col}_mean")
    ) %>% 
    ungroup()
  
}

# longitudinal clustering with mlm and k-means -------------------------------

longclus_mlm_kmeans_trad <- function(df, x,
                                     numclus, numruns,
                                     model,
                                     # model.diagcov=TRUE,
                                     trends='gcm',
                                     seed){
  
  
  
  df <- df %>% arrange(id)
  set.seed(seed)
  ids <- unique(df$id)
  #all_x <- df %>% select({{x}}) %>% as_vector()
  #x_time <- unique(all_x)
  #dt_time <- data.table(time=x_time)
  
  gcm <- eval(substitute(lmer(formula = model, data = df, REML = FALSE)))
  
  # need to substitute model.fixed and model.random otherwise predictY throws errors
  #gcm <- eval(substitute(hlme(fixed=model.fixed, random=model.random, subject='id', ng=1, idiag=model.diagcov, data=df)))
  Pred <- predict(gcm, df)
  # stopifnot(names(Pred) == as.integer(df$id)) #ensure proper order assumption for later code
  R <- ranef(gcm)$id %>% scale
  
  #capture.output(gcmsum <- summary(gcm)) %>% invisible #suppress printouts
  coef_fe <- t(fixef(gcm)) %>% as.data.table %>% #.[1,] %>%
    setnames(paste0('fixed_', names(.))) 
  names(coef_fe)[1] <- str_remove_all(names(coef_fe)[1],"\\(|\\)")
  coef_re <- data.table(id=ids, ranef(gcm)$id)
  coefs_traj <- cbind(id=ids, coef_fe[rep(1, nrow(coef_re)), ], coef_re[,-'id'])
  pred_traj <- data.table(df %>% select(id, {{x}}), Pred)
  
  #cld <- clusterLongData(traj=R, idAll=ids, time=seq_len(ncol(R)), varNames='Value')
  #par <- parALGO(saveFreq=1e99, scale=TRUE, startingCond='kmeans++')
  km_mlm <- kmeans(R, centers = numclus, nstart = numruns)
  #kml(cld, nbClusters=numclus, nbRedrawing=numruns, toPlot='none', parAlgo=par)
  #modk <- slot(cld, paste0('c', numclus))[[1]]
  clusters <- km_mlm$cluster
  #coefs_traj[, Cluster := clusters]
  pred_traj[, Cluster := clusters[id]]
  #postProbs <- modk@postProba
  #colnames(postProbs) <- levels(clusters)
  # compute centers
  #centers <- coefs_traj[, lapply(.SD, mean), keyby=Cluster, .SDcols=setdiff(colnames(coefs_traj), c('id', 'Cluster'))]
  # compute trends
  dt_trends <- pred_traj %>% 
    group_by(Cluster, {{x}}) %>% 
    summarise(Value = mean(Pred), .groups = "drop")
  dt_trends$Cluster <- as.factor(dt_trends$Cluster)
  # dt_trends <- pred_traj[, .(Value=mean(Pred)), by=.(Cluster, x)]
  # store results
  result <- list(df, clusters = tibble(id = ids, cluster = as.factor(clusters)), 
                 centers = km_mlm$centers, 
                 dt_trends=dt_trends, #start=start, postProbs=postProbs, 
                 numclus=numclus,
                 wss = km_mlm$tot.withinss)
  #result$BIC = modk@criterionValues['BIC']
  #result$AIC = modk@criterionValues[['AIC']]
  #result$converged = gcm$conv == 1
  #result$model$centers = centers
  
  return(result)
}

## plotting longitudinal clusters --------------------------------------------

plot_clusters <- function(df, clust_obj, x, y, group_var){
  
  clust_means <- clust_obj[["dt_trends"]] %>% 
    as_tibble() %>% 
    rename(cluster = Cluster)
  
  df %>% 
    left_join(clust_obj$clusters,
              by = "id") %>% 
    ggplot(aes(x = {{x}}, y = {{y}}, group = {{group_var}})) +
    geom_line(aes(colour = cluster), alpha = 0.5) +
    geom_line(data = clust_means, aes(x = {{x}}, y = Value,
                                      group = cluster,
                                      colour = cluster),
              size = 3) +
    scale_colour_brewer(palette = "Dark2") 
  
}

plot_clusters_facet <- function(df, clust_obj, x, y, group_var){
  
  clust_means <- clust_obj[["dt_trends"]] %>% 
    as_tibble() %>% 
    rename(cluster = Cluster)
  
  df %>% 
    left_join(clust_obj$clusters,
              by = "id") %>% 
    ggplot(aes(x = {{x}}, y = {{y}}, group = {{group_var}})) +
    geom_line(alpha = 0.5, colour = "darkgrey") +
    geom_line(data = clust_means, aes(x = {{x}}, y = Value,
                                      group = cluster),
              size = 2) +
    facet_wrap(~cluster) 
  
}

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
