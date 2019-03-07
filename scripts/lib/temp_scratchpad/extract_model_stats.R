temp_mean =  lapply(  GRL_splitbymodelkmer , function(x) x[1]$model_mean ) 
temp_sd   =  lapply(  GRL_splitbymodelkmer , function(x) x[1]$model_stdv ) 

temp_stats <- cbind( unlist(temp_mean), unlist(temp_sd) )
colnames( temp_stats ) <- c( "mean", "std_dev")


current_mean_observed <- lapply( GRL_splitbymodelkmer, function(x) mean(x$event_mean) )
current_sd_observed   <- lapply( GRL_splitbymodelkmer, function(x) mean(x$event_stdv) )

current_observed_stats <-  cbind( unlist(current_mean_observed), unlist(current_sd_observed) )
colnames( current_observed_stats ) <- c( "mean", "std_dev")
