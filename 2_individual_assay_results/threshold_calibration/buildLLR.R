# Downloaed from https://github.com/jweile/maveLLR/blob/master/R/buildLLR.R
# Copyright (C) 2020  Jochen Weile, Roth Lab
#
# This file is part of maveLLR
#
# maveLLR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# maveLLR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with maveLLR  If not, see <https://www.gnu.org/licenses/>.

#' Build LLR function using kernel density estimation
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernel <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=0.0001) {
  
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  
  #generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(posScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  #measure the hypothetical density of an outlier point
  minDensPos <- kdensity(c(posScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #do the same for the negative distribution
  pseudo <- min(negScores) - 10*ifelse(is.numeric(bw),bw,0.1)
  minDensNeg <- kdensity(c(negScores[-1],pseudo),bw=bw,kernel=kernel)(pseudo)
  #whichever outlier density is higher will serve as our uniform prior
  minDens <- outlierSuppression * max(minDensPos,minDensNeg)
  
  llrFun <- function(score) sapply(score,function(s)
    # log10( max(posDens(s),minDens) / max(negDens(s),minDens) )
    log10(max(posDens(s),minDens)) - log10(max(negDens(s),minDens))
  )
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}

#' Build LLR function using kernel density estimation - experimental version
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param bw the bandwith to use for kernel density estimation (see kdensity package for more info).
#'   This can either be a numerical value or the name of the algorithm used to automatically choose one.
#' @param kernel the type of kernel to use (see kdensity package for more info)
#' @param 
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.kernel(posScores,negScores,bw=0.1,kernel="gaussian")
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.kernelExperimental <- function(posScores, negScores, bw=0.1, kernel="gaussian", outlierSuppression=1) {
  ## modify the original outlierSuppression=0.0001 to outlierSuppression=1
  library(kdensity)
  
  posDens <- kdensity(posScores,bw=bw,kernel=kernel)
  negDens <- kdensity(negScores,bw=bw,kernel=kernel)
  refDens <- kdensity(c(posScores,negScores), bw=bw, kernel=kernel)
  
  # generate a dummy outlier point far away from the rest of the distribution
  # so we can use it to measure it's hypothetical density
  pseudo <- min(c(posScores, negScores)) - 10*ifelse(is.numeric(bw),bw,0.1)
  # measure the hypothetical density of an outlier point
  pseudoDens <- kdensity(c(posScores[-1],negScores,pseudo),bw=bw,kernel=kernel)(pseudo)
  priorWeight = pseudoDens * outlierSuppression
  
  llrFun <- function(score) sapply(score,function(s) {
    rawLLR = log10(posDens(s)) - log10(negDens(s))
    weight = refDens(s) / (refDens(s)+priorWeight)
    finalLLR = rawLLR * weight
  })
  
  return(list(llr=llrFun,posDens=posDens,negDens=negDens))
}


#' Build LLR function using Gaussian densities.
#'
#' @param posScores fitness scores in the positive reference set
#' @param negScores fitness scores in the negative reference set
#' @param spline logical; whether or not to apply spline monotonization. TRUE by default.
#'
#' @return a list containing the LLR function in log10 scale and the positive and negative density functions
#' @export
#'
#' @examples
#'
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
buildLLR.gauss <- function(posScores, negScores, spline=TRUE) {
  
  mpos <- mean(posScores)
  spos <- sd(posScores)
  mneg <- mean(negScores)
  sneg <- sd(negScores)
  
  posDens <- function(x) dnorm(x,mpos,spos)
  negDens <- function(x) dnorm(x,mneg,sneg)
  
  llrFun <- function(score) log10(posDens(score)/negDens(score))
  
  if (spline) {
    minPoint <- optimize(llrFun,interval=c(mpos, qnorm(0.999,mneg,sneg)),maximum=FALSE)
    if (minPoint$minimum < mneg) {
      minPoint$minimum <- Inf
    }
    maxPoint <- optimize(llrFun,interval=c(qnorm(0.001,mpos,spos),mneg),maximum=TRUE)
    if (maxPoint$maximum > mpos) {
      maxPoint$maximum <- -Inf
    }
    
    llrSpline <- function(scores) {
      sapply(scores,function(score) {
        if (is.na(score)) {
          NA
        } else if (score > minPoint$minimum) {
          minPoint$objective
        } else if (score < maxPoint$maximum) {
          maxPoint$objective
        } else {
          llrFun(score)
        }
      })
    }
    
    return(list(llr=llrSpline,posDens=posDens,negDens=negDens))
    
  } else {
    return(list(llr=llrFun,posDens=posDens,negDens=negDens))
  }
  
}


#' Draw a plot that summarizes the densities and LLR calculated by other functions
#'
#' @param scores the numerical scores in the Mave map
#' @param llrFun the LLR function
#' @param posDens the density function of the positive reference set
#' @param negDens the density function of the negative reference set
#' @param posScores the numerical scores in the postive reference set
#' @param negScores the numerical scores in the negative reference set
#'
#' @return nothing
#' @export
#'
#' @examples
#' posScores <- rnorm(50,0.2,0.3)
#' negScores <- rnorm(30,0.8,0.2)
#' llr <- buildLLR.gauss(posScores,negScores)
#' drawDensityLLR(c(posScores,negScores),llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#'
drawDensityLLR <- function(scores, llrFun, posDens, negDens, posScores, negScores, prior=0.1) {
  
  llrTs <- llrThresholds(optiLLR(prior))
  
  opar <- par(mfrow=c(2,1))
  xlim <- range(scores,na.rm=TRUE,finite=TRUE)
  ymax <- max(c(
    posDens(seq(xlim[[1]],xlim[[2]],length.out = 100)),
    negDens(seq(xlim[[1]],xlim[[2]],length.out = 100))
  ))
  par(mar=c(.1,4,1,1))
  xs <- seq(xlim[[1]],xlim[[2]],length.out=200)
  ys <- llrFun(xs)
  ylim <- range(ys,na.rm=TRUE,finite=TRUE)
  plot(NA,type="n",xlim=xlim,ylim=ylim,
       axes=FALSE,xlab="",ylab="LLR"
  )
  drawThresh <- function(t,col,label) {
    if (t >= 0) {
      if (t < ylim[[2]]) {
        rect(xlim[[1]],t,xlim[[2]],ylim[[2]],col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=3,cex=0.8)
      }
    } else {
      if (t > ylim[[1]]) {
        rect(xlim[[1]],ylim[[1]],xlim[[2]],t,col=col,border=NA)
        text(.1*xlim[[2]]+.9*xlim[[1]],t,label,pos=1,cex=0.8)
      }
    }
  }
  drawThresh(llrTs[["patho.support"]],"lemonchiffon","Patho. support.")
  drawThresh(llrTs[["patho.moderate"]],"lightgoldenrod1","Patho. moderate")
  drawThresh(llrTs[["patho.strong"]],"goldenrod1","Patho. strong")
  drawThresh(llrTs[["patho.vstrong"]],"indianred1","Patho. very str.")
  drawThresh(llrTs[["benign.support"]],"lemonchiffon","Benign support.")
  drawThresh(llrTs[["benign.strong"]],"goldenrod1","Benign strong")
  
  lines(xs,ys)
  # plot(llrFun,from=xlim[[1]],to=xlim[[2]],xlim=xlim,axes=FALSE,xlab="",ylab="LLR")
  abline(h=0,col="gray",lty="dashed")
  axis(2)
  par(mar=c(5,4,.1,1))
  hist(scores,col="gray90",border=NA,freq=FALSE,main="",xlim=xlim,ylim=c(0,ymax),breaks=50)
  plot.function(posDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="firebrick3",lwd=2)
  plot.function(negDens,from=xlim[[1]],to=xlim[[2]],add=TRUE,col="darkolivegreen3",lwd=2)
  abline(v=posScores,col="firebrick3")
  abline(v=negScores,col="darkolivegreen3")
  par(opar)
  
  return(invisible(NULL))
}

optiLLR <- function(prior=0.1,posterior=0.9) {
  log10(posterior*(1-prior)/(prior*(1-posterior)))*4/3
}

llrThresholds <- function(LLRpvst=optiLLR(0.1),X=2) {
  c(
    patho.vstrong=LLRpvst,
    patho.strong=LLRpvst/X,
    patho.moderate=LLRpvst/(X^2),
    patho.support=LLRpvst/(X^3),
    benign.support=-LLRpvst/(X^3),
    benign.strong=-LLRpvst/(X^1)
  )
}


#
# mthfr <- read.csv("~/projects/mthfr/folate_response_model5.csv")
# variants <- mthfr[mthfr$type=="substitution","hgvs"]
# scores <- mthfr[mthfr$type=="substitution","m25.score"]
#
# idb <- 337
# prsnrs <- read.csv("~/projects/mthfr/mthfr_prsNrs.csv")
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start < idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start < idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores,bw=0.1,kernel="gaussian")
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
#
# posref <- prsnrs[prsnrs$reference=="positive" & prsnrs$start > idb,"hgvs"]
# negref <- prsnrs[prsnrs$reference=="negative" & prsnrs$start > idb,"hgvs"]
# posScores <- na.omit(scores[variants %in% posref])
# negScores <- na.omit(scores[variants %in% negref])
# llr <- buildLLR.kernel(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)
# llr <- buildLLR.gauss(posScores, negScores)
# drawDensityLLR(scores,llr$llr,llr$posDens,llr$negDens,posScores,negScores)