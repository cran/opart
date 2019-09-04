## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(microbenchmark)
library(data.table)
library(directlabels)
library(ggplot2)

## ------------------------------------------------------------------------
data(neuroblastoma, package="neuroblastoma")
selectedData <- subset(neuroblastoma$profiles, profile.id=="1" & chromosome=="1")
result <- opart::opart_gaussian(selectedData$logratio, 1)
str(result)

## ----message=FALSE, warning=FALSE----------------------------------------
n.profiles <- 10
selected_ids <- unique(neuroblastoma$profiles$profile.id)
profile_ids_vec <- head(selected_ids, 10)
selected_data <- neuroblastoma$profiles

#initializing lists for storing the results of 3 algorithms
fpop_res <- list()
opart_res <- list()
cpt_res <- list()

#counting the iterations
iters <- 0

for(profile_id in profile_ids_vec){
  current_data = subset(selected_data, profile.id == profile_id)
  
  fpop_res <- fpop::Fpop(current_data$logratio, iters + 1)
  
  cpt_res <- changepoint::cpt.mean(data = current_data$logratio, 
                        penalty = "Manual", 
                        method = "PELT", 
                        pen.value=iters + 1)
  
  opart_res <- opart::opart_gaussian(current_data$logratio, iters + 1)
  
  iters <- iters + 1
  
  if(!isTRUE(all.equal(fpop_res$t.est, opart_res$end.vec)) || 
     !isTRUE(all.equal(cpt_res@cpts, opart_res$end.vec)))
    break
}

## ------------------------------------------------------------------------
iters

## ------------------------------------------------------------------------
timing_list <- list() #for storing the results of microbenchmark

data.counts <- data.table(neuroblastoma$profiles)[, .(count=.N), by=.(profile.id, chromosome)]
data.counts <- data.counts[order(count)][as.integer(seq(1, .N, l=100))]

for(i in 1:nrow(data.counts)){
  current_data = subset(neuroblastoma$profiles, profile.id == data.counts[i, profile.id] 
                        & chromosome == data.counts[i, chromosome] )
  size <- length(current_data$logratio)
  
  timing <- microbenchmark(
  "fpop"={
      fpop::Fpop(current_data$logratio, 1)
    },
  "changepoint"={
      changepoint::cpt.mean(data = current_data$logratio, 
                            penalty = "Manual", 
                            method = "PELT", 
                            pen.value=1)
    },
   "opart"={
        opart::opart_gaussian(current_data$logratio, 1)
    }, times=5)
  
  timing_list[[paste(i)]] <- data.table(size, timing)
}

timing.dt <- do.call(rbind, timing_list)

## ------------------------------------------------------------------------
timing.dt$time <- timing.dt$time * (10^(-9))

## ----message = FALSE, fig.height = 5, fig.width = 10---------------------

lab.df <- data.frame(
seconds <- c(1,60,3600),
labels <- c("1 second", "1 minute", "1 hour"))

gg_runtime <- ggplot(data = timing.dt, aes(x = (size),y = (time),
                                  col = expr))+
geom_point()+
geom_smooth()+
  
geom_hline(aes(yintercept=(seconds)),
             data=lab.df,
             color="grey")+
geom_text(aes(100, (seconds), label=labels),
            data=lab.df,
            size=3,
            color="black",
            vjust=-0.5)+
scale_x_log10("size of data vector") + scale_y_log10("time(s)") 

direct.label(gg_runtime, "angled.boxes")

## ------------------------------------------------------------------------
timing_list <- list() #for storing the results of microbenchmark
selected_data <- rnorm(10000, mean=100)
opfit <- c()
fpfit <- c()
cpfit <- c()
for(penalty in 10^(-10:10)){
  timing <- microbenchmark(
    "opart"={
        opfit <- opart::opart_gaussian(selected_data, penalty)
    }, 
    "fpop"={
        fpfit <- fpop::Fpop(selected_data, penalty)
    },
    "changepoint"={
        cpfit <- changepoint::cpt.mean(data = selected_data, penalty = "Manual", 
               method = "PELT", pen.value=penalty)
    },times=5)
   
  timing$cpts <- 1:nrow(timing)
  timing[timing$expr=="opart",]$cpts <- length(opfit$end.vec)
  timing[timing$expr=="fpop",]$cpts <- length(fpfit$t.est)
  timing[timing$expr=="changepoint",]$cpts <- length(cpfit@cpts)
  timing_list[[paste(penalty)]] <- data.table(penalty, timing)
}
penalty_timing <- do.call(rbind, timing_list)

## ------------------------------------------------------------------------
penalty_timing$time <- penalty_timing$time * (10^(-9))

## ----message = FALSE, fig.height = 5, fig.width = 10---------------------
gg_penalty <- ggplot(data = penalty_timing, aes(x = penalty,y = time,
                                  col = expr, label=cpts))+
geom_point()+
geom_smooth()+
geom_label()+
scale_x_log10("penalty") + scale_y_log10("time(s)")+
  scale_colour_discrete(name  ="Package",
                        labels=c("opart", "fpop", "changepoint"))
gg_penalty

