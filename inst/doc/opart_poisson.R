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
sample_data <- rpois(n = 100, lambda = 10)
opfit <- opart::opart_poisson(data = sample_data, penalty = 1)
str(opfit)

## ----message=FALSE, warning=FALSE----------------------------------------
sfit <- Segmentor3IsBack::Segmentor(c(1,2,3,4), Kmax = 4)
sfit

## ------------------------------------------------------------------------
str(sfit)

## ------------------------------------------------------------------------
sim.seg <- function(seg.mean, size.mean=15){
seg.size <- rpois(1, size.mean)
rpois(seg.size, seg.mean)
}
set.seed(1)
seg.mean.vec <- c(1.5, 3.5, 0.5, 4.5, 2.5)
z.list <- lapply(seg.mean.vec, sim.seg)
(z.rep.vec <- unlist(z.list))

## ------------------------------------------------------------------------
count.df <- data.frame(
chrom="chrUnknown",
chromStart=0:(length(z.rep.vec)-1),
chromEnd=1:length(z.rep.vec),
count=z.rep.vec)

gg.count <- ggplot()+
xlab("position")+
geom_point(aes(
chromEnd, count),
shape=1,
data=count.df)
gg.count

## ------------------------------------------------------------------------
n.segs <- length(seg.mean.vec)
seg.size.vec <- sapply(z.list, length)
seg.end.vec <- cumsum(seg.size.vec)
change.vec <- seg.end.vec[-n.segs]+0.5
change.df <- data.frame(
changepoint=change.vec)
gg.change <- gg.count+
geom_vline(aes(
xintercept=changepoint, color=model),
data=data.frame(change.df, model="simulation"))+
scale_color_manual(
values=c(
simulation="black",
opart="green",
Segmentor3IsBack="blue"))
gg.change

## ------------------------------------------------------------------------
opfit <- opart::opart_poisson(z.rep.vec, 10.5)
opend.vec <- opfit$end.vec

chromStart <- c(1, head(opend.vec, -1) + 1)

opfit.segments <- data.frame(chromStart <- (chromStart), 
                       chromEnd <- (opend.vec))

names(opfit.segments) <- c("chromStart", "chromEnd")

chromMean <- c()

for (i in 1:nrow(opfit.segments)){
  s <- opfit.segments[i,]
  mu <- mean(z.rep.vec[s["chromStart"][1,] : s["chromEnd"][1,]])
  chromMean[paste(i)] <- mu
}

opfit.segments["chromMean"] <- chromMean

gg.change+
geom_segment(aes(
chromStart, chromMean, xend=chromEnd, yend=chromMean, color=model),
data=data.frame(opfit.segments, model="opart"))

## ------------------------------------------------------------------------
sfit <- Segmentor3IsBack::Segmentor(data=z.rep.vec, Kmax=length(opfit$end.vec))
sfitend.vec <- as.numeric(sfit@breaks[length(opfit$end.vec),])

chromStart <- c(1, head(sfitend.vec, -1) + 1)

sfit.segments <- data.frame(chromStart <- (chromStart), 
                       chromEnd <- (sfitend.vec))

names(sfit.segments) <- c("chromStart", "chromEnd")

chromMean <- c()

for (i in 1:nrow(sfit.segments)){
  s <- sfit.segments[i,]
  mu <- mean(z.rep.vec[s["chromStart"][1,] : s["chromEnd"][1,]])
  chromMean[paste(i)] <- mu
}

sfit.segments["chromMean"] <- chromMean

gg.change+
geom_segment(aes(
chromStart, chromMean, xend=chromEnd, yend=chromMean, color=model),
data=data.frame(sfit.segments, model="Segmentor3IsBack"))

