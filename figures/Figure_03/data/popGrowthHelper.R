
# Helper functions for estimating cell growth and fitness along time course
# =========================================================================

##############################
#        Eran Kotler         #
#  Last updated: 2020-10-30  #
##############################

library(ggplot2)
library(tidyverse)
library(gridExtra)

growthCurve <- function(d, loess.span=3, deg.poly=2, deg.fr=5, title_pfx="", fit.type="loess", x.to.pred=NA, out_f=NA){
  # Plot growth curve, and fit approximation (either "loess", "poly" or "splines"). Returns value or derivative at x.to.pred
  # d is a datafrme with Total.viable.cells and Cells.seeded for each NumDays (columns). Comes from Kasper's organoid logbooks
  d <- as.data.table(d)
  d <- d[!is.na(Total.viable.cells),] # original, works
    
  d[d$Total.viable.cells==0,"Total.viable.cells"] <- 1 # Avoid problems in LOESS in cultures where population drops to 0 (e.g. P1 LATE D3 R2)

    
  setorder(d, NumDays)
  d[ , FCgrowth := Total.viable.cells / shift(Cells.seeded, 1L, type="lag")] # Growth (FC in size) from previous time point (point-wise increase in pop size)
  
  # Calc cumulative population size as [pop size at previous time point] X [FC growth from orevious time point to current] 
  d[1, "Cum.pop.size"] = d[1,"Total.viable.cells"]
  for (t in seq(2,nrow(d))){
    prev.pop <- d[t-1, "Cum.pop.size"]
    curr.FCgrowth <- d[t, "FCgrowth"]
    d[t, "Cum.pop.size"] <- curr.FCgrowth * prev.pop
  }
  d$Log.pop.size <- log10(d$Cum.pop.size)
  
  # fit splines and get derivative:
  if (fit.type=="loess"){
    if (!is.na(x.to.pred)){ # print out derivative at requested point
      s <- sprintf("Slope at %i: %f",x.to.pred,
                   loess.approximation(d$NumDays, d$Log.pop.size, loess.span=loess.span, 
                                       show.plot=F, x.to.pred=x.to.pred, return_deriv=T))
      print(s)}
    p <- loess.approximation(d$NumDays, d$Log.pop.size, loess.span=loess.span, show.plot=T, title_pfx=title_pfx)
  }else{if (fit.type=="poly"){
    p <- polynom.approximation(d$NumDays, d$Log.pop.size, deg.poly=deg.poly, show.plot=T, title_pfx=title_pfx)
  }else{if (fit.type=="splines"){   
    p <- splines.approximation(d$NumDays, d$Log.pop.size, deg.fr=deg.fr, show.plot=T, title_pfx=title_pfx)    
  }}}
    
  if(!is.na(out_f)){ # Save plots to file
      ggsave(file=out_f, p, width = 9, height = 9)
      print(paste("Saved plot to:", out_f))
    }
  return(p)
}

splines.approximation <- function(x, y, x.to.pred=0, return_deriv=FALSE, deg.fr=NA, show.plot=FALSE, title_pfx=""){
  # Fit smooth spline to growth curve and return approximation (value or 1st derivative) at requested point 
  # x.to.pred: x value for which to return the approximated y according to the fit
  # return_deriv: If TRUE returns 1st derivative instead of estimated value at point
  # deg.fr: degrees of freedom of smooth splines (if NA determined automatically)
  # If show.plot=FALSE returns a predicted value for x.to.pred)
  
  xx <- unique(sort(c(seq(0, 1000, by = 1), unique(x)))) # grid for x axis
  spl <- smooth.spline(x, y, df = deg.fr)
  # spl <- smooth.spline(x, y, cv = TRUE) # choose df by CV (tends to over-fit here)
  der <- ifelse(return_deriv==TRUE, yes = 1, no=0) # derivative to return
  pred <- predict(spl, x.to.pred, deriv = der)
  if (show.plot){     
    
    title.str <- paste(title_pfx, sprintf("smooth spline with df = %f",spl$df))
    p1 <-ggplot(data=tibble(x=x, y=y, splines=predict(spl,x)[[2]])) + #plot data and spline fit
      geom_point(aes(x=x, y=y), pch=20, color="blue", size=2) +
      geom_line(aes(x=x, y=splines) ) +
      ggtitle(title.str)+
      ylab("Estimated population size [log10(#cells)]") + xlab("Time [days]") +
      ylim(0,60) + xlim(100,600)
    
    p2 <-ggplot(data=tibble(x=x, y=predict(spl, deriv = 1)[[2]])) + # plot derivs
      geom_point(aes(x=x, y=y), pch=20, color="red", size=2) +
      ylab("Derivative") + xlab("Time [days]")+
      ylim(-0.02,0.15) + xlim(100,600)
    
    g <- grid.arrange(p1, p2, nrow=2)
    return(g)
  }
  return(pred$y) # If not plotting - return 1st derivative at point requested
}


eval.deriv <- function(fit.coefs, x.point){
  # Calc deriv of 2nd degree polynomial at given point
  if (length(fit.coefs)>3){
    stop("NOT IMPLEMENTED FOR POLYNOMIALS OF HIGHER DEGREE THAN 2")
  }
  der <- fit.coefs[2] + (2*x.point*fit.coefs[3])
  as.numeric(der)
}


polynom.approximation <- function(x, y, x.to.pred=0, return_deriv=FALSE, deg.poly=2, show.plot=FALSE, title_pfx=""){
  # Fit polynomial to growth curve and return approximation (value or 1st derivative) at requested point 
  # x.to.pred: x value for which to return the approximated y according to the fit
  # return_deriv: If TRUE returns 1st derivative instead of estimated value at point
  # deg.fr: degrees of freedom of smooth splines (if NA determined automatically)
  # If show.plot=FALSE returns a predicted value for x.to.pred)
  fit <- lm(y~poly(x=x, degree=deg.poly, raw=TRUE)) # fit polynomial
  xlims <- range(x)
  out <- ifelse(return_deriv==TRUE, yes = eval.deriv(fit$coefficients, x.to.pred), predict(fit, newdata=data.frame(x=x.to.pred)))
  
  if (show.plot){
    new <- data.frame(x=seq(xlims[1], xlims[2])) # x values over which to plot fit
    preds <- predict(fit, newdata=new) # predictions (approximations)
    derivs <- sapply(x, function(x.to.pred) eval.deriv(fit$coefficients, x.to.pred))
    title.str <- paste(title_pfx, "2nd deg. polynimoial")

    p1 <-ggplot(data=tibble(x=x, y=y)) + #plot data and spline fit
      geom_point(aes(x=x, y=y), pch=20, color="blue", size=2) +
      geom_line(data=tibble(x=new$x, y=preds), aes(x=x, y=y) ) +
      ggtitle(title.str)+
      ylab("Estimated population size [log10(#cells)]") + xlab("Time [days]") #+

    p2 <-ggplot(data=tibble(x=x, y=derivs)) + # plot derivs
      geom_point(aes(x=x, y=y), pch=20, color="red", size=2) +
      ylab("Derivative") + xlab("Time [days]") #+
      # ylim(-0.02,0.15) + xlim(100,600)
    g <- grid.arrange(p1, p2, nrow=2)
    return(g)
  }
  return(out) # If not plotting - return 1st derivative at point requested
}


loess.approximation <- function(x, y, x.to.pred=NA, return_deriv=FALSE, loess.span=1, show.plot=FALSE, title_pfx=""){
  # Fit local regression to growth curve and return approximation (value or 1st derivative) at requested point 
  # x.to.pred: x value for which to return the approximated y according to the fit
  # return_deriv: If TRUE returns numerically approximaterd 1st derivative instead of estimated value at point
  # loess.span: a parameter for LOESS (controls degree of smoothing - low value-> more local)
  # If show.plot=FALSE returns a predicted value for x.to.pred)
  
  grid.res <- 0.01 # resolution of prediction grid (and of slope)
  x.grid <- seq(range(x)[1]-grid.res, range(x)[2]+grid.res, grid.res)
  l <- loess(y~x, span=loess.span) # fit loess
  preds <- predict(l, data.frame(x=x.grid))
  
  # Numerical estimation of the slope:
  sl <- c(NA, diff(preds,lag=2)/(x.grid[3]-x.grid[1]), NA) # dy/dx
  
  if(return_deriv==TRUE){ # If deriv is not defined at exact point - take the next point on grid (for edges)
    exact <- sl[x.grid==x.to.pred]
    out <- ifelse(!is.na(exact), 
                  yes=exact, 
                  no=ifelse(!is.na(sl[x.grid==x.to.pred+grid.res]), 
                            yes = sl[x.grid==x.to.pred+grid.res],
                            no = sl[x.grid==x.to.pred-grid.res]))
  } else { out <- predict(l, newdata=data.frame(x=x.to.pred)) }
  
  if (show.plot){
    title.str <- paste(title_pfx, "LOESS")
    x_lims = c(0, max(x)+20)
    p1 <-ggplot(data=tibble(x=x, y=y)) + #plot data and spline fit
      geom_point(aes(x=x, y=y), pch=20, color="blue", size=2) +
      geom_line(data=tibble(x=x.grid, y=preds), aes(x=x, y=y) ) +
      
      ggtitle(title.str)+
      ylab("Estimated population size [log10(#cells)]") + xlab("Time [days]") +
      coord_cartesian(ylim=c(0,100), xlim=x_lims) + # c(100,600)
    if (!is.na(x.to.pred)){p1 <- p1+geom_vline(aes(xintercept=x.to.pred), linetype = "dashed")}
    
    p2 <-ggplot(data=tibble(x=x.grid, y=sl)) + # plot derivs
      geom_line(aes(x=x, y=y), pch=20, color="red") +
      ylab("Derivative") + xlab("Time [days]") +
      coord_cartesian(ylim=c(-0.02,0.15), xlim=x_lims) + # c(100,600)
    if (!is.na(x.to.pred)){p2 <- p2+geom_vline(aes(xintercept=x.to.pred), linetype = "dashed")}
    
    # remove grey background & grid:
    p1 <- p1 + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) # theme_bw()
    p2 <- p2 + theme(panel.background = element_blank()) + theme(axis.line = element_line(colour = "black")) # theme_bw()

    g <- grid.arrange(p1, p2, nrow=2)
    return(g)
  }
  return(as.numeric(out)) # If not plotting - return 1st derivative at point requested
}

