reconstruct <-
function(WT, my.series = 1, lvl = 0, 
      only.coi = FALSE, 
      only.sig = TRUE, siglvl = 0.05, 
      only.ridge = FALSE, 
      sel.period = NULL, sel.lower = NULL, sel.upper = NULL,  
      rescale = TRUE,
      plot.waves = FALSE, plot.rec = TRUE, 
      lty = 1, lwd = 1, col = 1:2, ylim = NULL,
      show.legend = TRUE, 
      legend.coords = "topleft", legend.horiz = FALSE, legend.text = NULL,
      label.time.axis = TRUE, 
      show.date = FALSE, date.format = NULL, date.tz = NULL,
      timelab = NULL, timetck = 0.02, timetcl = 0.5,
      spec.time.axis = list(at = NULL, labels = TRUE, 
                            las = 1, hadj = NA, padj = NA),
      main.waves = NULL, main.rec = NULL, main = NULL, 
      lwd.axis = 1, 
      verbose = TRUE) {
         
    if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    } 
    
#     lwd.axis = 0.25
    
    series.data = WT$series
    
    ####################################
    ## Identify the scenario
    ####################################
    
    if (class(WT) == 'analyze.wavelet') {
    
        out("Your input object class is 'analyze.wavelet'...\n") 
        
        my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[2], names(series.data)[1]) 
        
        Wave  = WT$Wave
        Power = WT$Power
        Power.pval = WT$Power.pval
        Ridge = WT$Ridge     
        
    } 
    if (class(WT) == 'analyze.coherency') {   
    
        out("Your input object class is 'analyze.coherency'...\n")   
        
        if (is.numeric(my.series)) { 
            if (!is.element(my.series,c(1,2))) { stop("Please choose either series number 1 or 2!") }
            my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[my.series+1], names(series.data)[my.series])  
        }
        
        ind = which( names(series.data) == my.series ) 
        which.series.num = ifelse(names(series.data)[1] == 'date', ind-1, ind)
        if (!is.element(which.series.num, c(1,2))) { stop("Your series name is not available, please check!") }
      
        if (which.series.num == 1) {
            Wave  = WT$Wave.x
            Power = WT$Power.x
            Power.pval = WT$Power.x.pval
            Ridge = WT$Ridge.x
        }
        if (which.series.num == 2) {
            Wave  = WT$Wave.y
            Power = WT$Power.y
            Power.pval = WT$Power.y.pval
            Ridge = WT$Ridge.y
        }      
        
    }   
           
    out(paste("Your time series '", my.series, "' will be reconstructed...", sep=''), '\n')   
   

    ####################################
    ## Prepare reconstruction components
    ####################################
    
    out("Starting the reconstruction process...\n")
    
    nc = WT$nc    
    nr = WT$nr

    dt = WT$dt
    dj = WT$dj
    
    Scale = WT$Scale
    Period = WT$Period   
       
    loess.span = WT$loess.span
   
   
    rec.waves = matrix(0, nrow=nr, ncol=nc)
    for (s.ind in seq_len(nr)) {
         rec.waves[s.ind,] = (Re(Wave[s.ind,])/sqrt(Scale[s.ind]))*dj*sqrt(dt)/(pi^(-1/4)*0.776)
    }
    
    # select minimum level?
    rec.waves = rec.waves * (Power >= lvl)
    comment.lvl = paste('minimum power level: ',lvl, sep='')
    
    # select ridge?
    if (only.ridge == T) {      
        rec.waves = rec.waves * Ridge  
        rec.waves[Ridge == 0] = NA
    }  
    comment.ridge = paste('only ridge: ', only.ridge, sep='')
    

    # use significant parts only?
    if (only.sig == T) {
        if (!is.null(Power.pval)) { 
            rec.waves = rec.waves * (Power.pval < siglvl)
            rec.waves[Power.pval >= siglvl] = NA
        }      
    }  
    if (only.sig == F) { siglvl = NA }
    comment.sig = paste('significance level: ',siglvl,sep='')
    
    # only cone of influence for reconstruction? 
    if (only.coi == T) {
        for (i in (1:nc)) {
             for (s.ind in seq_len(nr)) {
                  if (Scale[s.ind] > 2^WT$coi.2[i]) {rec.waves[s.ind, i] = NA}
             }
        } 
    }  
    comment.coi = paste('only coi: ', only.coi, sep='') 
    
    # so far
    rnum.used = which(rowSums(rec.waves, na.rm=T)!=0)
    comment.periods = 'period: all relevant'
    
       
    # sel.period available?
    if (length(sel.period) != 0) {
    
        sel.rnum = numeric()
    
        # nearest available period:
        for (i in (1:length(sel.period))) {
             sel.rnum = union(sel.rnum, which(abs(Period-sel.period[i]) == min(abs(Period-sel.period[i])))) 
        }  
        
        rec.waves = rec.waves[sel.rnum,]
        comment.periods = paste('period: ',paste(as.character(round(Period[sel.rnum], 1)), collapse=', '), sep='')
           
        if (length(sel.rnum) == 1) {
            rec.waves = t(rec.waves)
        } 
        
        rnum.used = intersect(rnum.used, sel.rnum)
    } 
    
    # in case, sel.period is not available, refer to sel.upper, sel.lower
    if (length(sel.period) == 0 & ((length(sel.lower) != 0) | (length(sel.upper) != 0))) {   
    
       # in case, sel.lower (sel.upper) is not available, use minimum (maximum) period 
       if (length(sel.lower) == 0) {sel.lower = min(Period)}
       if (length(sel.upper) == 0) {sel.upper = max(Period)}             
       # in case, sel.lower > sel.upper, interchange
       if (sel.lower > sel.upper) { 
           sel.lower.h = sel.lower
           sel.lower = sel.upper
           sel.upper = sel.lower.h
       }
       sel.rnum = which(((Period >= sel.lower) & (Period <= sel.upper))) 
       rec.waves = rec.waves[sel.rnum,]
       
       # selected band / range of periods   
       sel.period.band  = Period[sel.rnum]  
       sel.period.range = as.character(round(sel.period.band, 1))
       if (length(sel.rnum) > 1) {
           sel.period.range = paste(round(range(sel.period.band), 1), collapse=' - ')
       }   
       if (length(sel.rnum) == 1) {rec.waves = t(rec.waves)}
       comment.periods = paste('period: ', sel.period.range, sep='') 
       
       rnum.used = intersect(rnum.used, sel.rnum)
       
    }
         
        
    
    ####################################
    ## Compute reconstructed series
    ####################################   
       
    # reconstructed time series
    x.r  = colSums(rec.waves, na.rm=T)
    
    # retrieve original time series
    x    = series.data[[my.series]]
    
    # rescale the reconstructed time series?
    if (rescale == T) {
        x.r  = (x.r-mean(x.r))*sd(x)/sd(x.r) + mean(x) 
    }   
    
    
    ####################################
    ## Plottings check: time axis default or not default?
    ####################################
    
    # date parameters
  
    if (is.null(date.format)) { date.format = WT$date.format }
    if (is.null(date.tz)) { date.tz = ifelse(is.null(WT$date.tz),"",WT$date.tz) }
  
    # individual time axis parameters
  
    if (!is.list(spec.time.axis))       spec.time.axis = list()
    if (is.null(spec.time.axis$at))     spec.time.axis$at = NULL
    if (is.null(spec.time.axis$labels)) spec.time.axis$labels = T
    if (is.null(spec.time.axis$las))    spec.time.axis$las = 1
    if (is.null(spec.time.axis$hadj))   spec.time.axis$hadj = NA
    if (is.null(spec.time.axis$padj))   spec.time.axis$padj = NA
  
    # initialize warning indicators for the case of individual time axis specification
    # warning: reset to time axis default?
    time.axis.warning = F 
    # warning: NAs among time axis tick marks?
    time.axis.warning.na = F
    # warning: calendar dates not chronological and/or format not standard unambiguous 
    chronology.warning = F
    
    if ( (!is.null(spec.time.axis$at)) & (label.time.axis==F) )  warning("\nPlease set label.time.axis = TRUE to make time axis specification effective.", immediate. = TRUE)
   
    ##############################
    
    # label time axis ? 
    if (label.time.axis == T) {
  
      # check if there is a user-defined time axis specification
      time.axis.default = ( is.null(spec.time.axis$at) )
  
      # check chronology and format of calendar date in case show.date = T
      if (show.date==T) { 
      
           if (is.element('date',names(series.data))) { my.date = series.data$date } else { my.date = rownames(series.data) }
           
           if (is.null(date.format)) { 
               
                   chronology.warning = inherits(try(as.Date(my.date, tz=date.tz), silent=T),'try-error')
                   
                   if (!chronology.warning) {
                       my.date = as.Date(my.date, tz=date.tz) 
                       chronology.warning = ifelse( sum(is.na(my.date))> 0, TRUE, sum(diff(my.date, tz=date.tz)<0)>0 )
                   }
                   
           }
           if (!is.null(date.format)) {
               
                   chronology.warning = inherits(try(as.POSIXct(my.date, format=date.format, tz=date.tz), silent=T),'try-error')
               
                   if (!chronology.warning) {
                       my.date = as.POSIXct(my.date, format=date.format, tz=date.tz) 
                       chronology.warning = ifelse( sum(is.na(my.date))> 0, TRUE, sum(diff(my.date, tz=date.tz)<0)>0 )
                   }
            }
            if (chronology.warning) { 
                show.date = F
                time.axis.default = T
                timelab = 'index'
            }
            
      } # end if (show.date==T)      
       
      # first: check conditions of time axis specification
      # time.axis.warning = TRUE leads to default reset
      
      if ( !time.axis.default ) { 
      
           # 1. According to show.date: Are tick marks (spec.time.axis$at) interpretable as numbers or dates?
           # show.date = F requires numeric values
           # show.date = T requires dates
           # mismatch results in default reset
           
           # case show.date==FALSE
           # Are tick marks (spec.time.axis$at) interpretable as numbers?
                      
           # is.numeric returns TRUE iff there is no character value but at least one non-NA value which is interpretable as number
           # (further logical entries T and F are interpreted as 1 and 0 which could be accepted here, non-positive values as well,
           # since values are matched with the true range of values in the axis command)
           
           if (show.date==F) { time.axis.warning = (!is.numeric(spec.time.axis$at)) }
           
           # case show.date==TRUE
           # Are tick marks (spec.time.axis$at) interpretable as dates using given date format and time zone specifications?
           
           # checking the first value given in spec.time.axis$at:
           # inherits is TRUE (and thus gives a time axis warning) iff the first value cannot be formatted, regardless of subsequent values
           # date formatting may produce NAs, however, therefore checking is necessary whether there are NAs ONLY
           # (NA tick marks are omitted automatically by the axis command as long as there are non-NAs)
           
           if (show.date==T) { 
           
               if (is.null(date.format)) { 
                   time.axis.warning = inherits(try(as.Date(spec.time.axis$at, tz=date.tz), silent=T), 'try-error') 
                   if (!time.axis.warning) { time.axis.warning = (sum(!is.na(as.Date(spec.time.axis$at, tz=date.tz)))==0) }
               }
            
               if (!is.null(date.format)) { 
                   time.axis.warning = inherits(try(as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz), silent=T),'try-error') 
                   if (!time.axis.warning) { time.axis.warning =  (sum(!is.na(as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz)))==0) }
               }            
           
           } # end if (show.date==T)
           
           if (!time.axis.warning) { 
           # 2. Are labels (spec.time.axis$labels) appropriate? Either logical and single-valued or non-logical and of the same length as tick marks?
           # is.logical returns TRUE iff there is not any character or numeric value, but there could be NA, even if it is the only value.
           
             # if logical:  single-valued and non-NA? 
             if ( is.logical(spec.time.axis$labels) ) { 
                  time.axis.warning = ( ifelse(length(spec.time.axis$labels)!=1,TRUE, is.na(spec.time.axis$labels)) )
             }
           
             # if non-logical: do vectors of tick marks and labels have equal length?
             if (!is.logical(spec.time.axis$labels)) {
                 time.axis.warning = ( length(spec.time.axis$labels) != length(spec.time.axis$at) ) 
             }
           } # end if (!time.axis.warning)
           
      } # end if ( !time.axis.default )
       
      # default reset in case of a warning
      time.axis.default = (time.axis.default | time.axis.warning)
      
      if ( (is.null(timelab) & time.axis.default) | (!is.null(timelab) & time.axis.warning) ) {timelab=ifelse(show.date, 'calendar date','index')}
           
    } # end if (label.time.axis == T)
    
    ####################################
    ## Plottings
    ####################################
    
    ## plot of reconstruction waves
        
    if (plot.waves == T) {
    
        out("Reconstruction waves are being plotted...\n")
        
        if (is.null(main) == F) {main.waves = main}
        
        range.rec.waves = range(rec.waves, na.rm =T)
#         matplot(WT$axis.1, t(rec.waves), type = 'l', ylim = range.rec.waves,
        matplot(1:nc, t(rec.waves), type = 'l', ylim = range.rec.waves,
                main = main.waves, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
                xaxs = 'i', xaxt = 'n', 
                xlab = '', ylab = '')

       # label time axis ?   
       if (label.time.axis == T) {
           
           #####################
      
           if (show.date == F) {
           
               if (time.axis.default) {
               
                   A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = timetck, tcl = timetcl)
                   mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
                   
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   time.tick = spec.time.axis$at
                   time.tick.label = spec.time.axis$labels
                   
                   # NAs among tick marks?
                   which.na = which(is.na(time.tick))
                   # NA warning
                   if ( length(which.na)>0 ) { time.axis.warning.na = T } 
                   
                   # in case of NAs among tick marks:
                   # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
                   axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
                        las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                       
               } # end if (!time.axis.default)
               
               mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)
               
           } # end if (show.date == F)
      
           #####################
           
           if (show.date == T) {  
           
               par(new=TRUE)
               
               if (time.axis.default) {
               
                   # empty plot, but calendar axis (w/ plotting)
                   plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],length.out=nc), 
                        ylim = range.rec.waves,
                        type="n", 
                        xaxs = "i", 
                        yaxt='n', 
                        xlab="", ylab="",
                        lwd = lwd.axis, tck=timetck, tcl=timetcl, 
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                        
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   # empty plot, but calendar axis (w/o plotting)
                   plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],length.out=nc), 
                        ylim = range.rec.waves,
                        type="n", 
                        xaxs = "i", 
                        xaxt='n', yaxt='n', 
                        xlab="", ylab="")
                                           
                   # user-defined calendar axis
                   if (is.null(date.format)) { time.tick = as.Date(spec.time.axis$at, tz=date.tz) }
                   if (!is.null(date.format)) { time.tick = as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz) }
                   
                   time.tick.label = spec.time.axis$labels
                   
                   # NAs among tick marks?
                   which.na = which(is.na(time.tick))
                   # NA warning
                   if ( length(which.na)>0 ) { time.axis.warning.na = T } 
                   
                   # in case of NAs among tick marks:
                   # NA tick marks (and corresponding labels) are omitted automatically by the axis command
                   
                   if ( is.logical(time.tick.label) ) {
                        # then time.tick.label has length 1 by criterion and IS TRUE OR FALSE
                        if (time.tick.label == T) {time.tick.label = time.tick}
                   }   
          
                   axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
                        las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                    
               } # end if (!time.axis.default)             
               
               mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)  
               
           }  # end if (show.date == T)   
           
       } # end if (label.time.axis == T)                  
                
       
    } #end if (plot.waves == T)
    
    x.rec = cbind(series.data, x.r=x.r)
    my.rec.series = paste(my.series,'.r',sep='')
    colnames(x.rec) = c(colnames(series.data),my.rec.series)
    rownames(x.rec) = rownames(series.data)
        
    if (plot.waves & plot.rec) { 
        par(ask = T)
    }
    
    
    
    ## plot of reconstructed series
    
    if (plot.rec == T) {
          
       out("Original (detrended) and reconstructed series are being plotted...\n")
       
       if (is.null(main) == F) { main.rec = main }
       
       range.rec = range(x.rec[,c(my.series,my.rec.series)], na.rm=T)
       if (is.null(ylim)) {ylim = range.rec} 
#        matplot(WT$axis.1, x.rec[,c(my.series,my.rec.series)], type = 'l', ylim = ylim, lty = lty, lwd = lwd, col = col,
       matplot(1:nc, x.rec[,c(my.series,my.rec.series)], type = 'l', ylim = ylim, lty = lty, lwd = lwd, col = col,
               main = main.rec, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
               xaxs = 'i', xaxt = 'n', 
               xlab = '', ylab = '')
               
       par(ask = F)
       
       if (show.legend == T) {      
       
           if (is.null(legend.text)) { 
               legend.text = c(paste('original', ifelse(loess.span!=0, paste(' (detrended, span: ', loess.span, ')', sep=''),''), sep=''), 'reconstructed')
               }
               
           legend(legend.coords, horiz=legend.horiz, legend = legend.text, lty = lty, lwd = lwd, col = col, text.font = par()$font.lab, cex=par()$cex.lab)        
       }
       
       
       # label time axis ?   
       if (label.time.axis == T) {
           
           #####################
           
           if (show.date == F) {
           
               if (time.axis.default) {
               
                   A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = timetck, tcl = timetcl)
                   mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
                   
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   time.tick = spec.time.axis$at
                   time.tick.label = spec.time.axis$labels
                   
                   # NAs among tick marks?
                   which.na = which(is.na(time.tick))
                   # NA warning
                   if ( length(which.na)>0 ) { time.axis.warning.na = T } 
                   
                   # in case of NAs among tick marks:
                   # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
                   axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
                        las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                                     
               } # end if (!time.axis.default)
               
               mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)
               
           } # end if (show.date == F)
      
           #####################
           
           if (show.date == T) {  
           
               par(new=TRUE)
               
               if (time.axis.default) {
               
                   # empty plot, but calendar axis (w/ plotting)
                   plot(my.date, seq(ylim[1], ylim[2], length.out=nc), 
                        ylim = ylim,
                        type="n", 
                        xaxs = "i", 
                        yaxt='n', 
                        xlab="", ylab="",
                        lwd = lwd.axis, tck=timetck, tcl=timetcl, 
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)                         
                        
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   # empty plot, but calendar axis (w/o plotting)
                   plot(my.date, seq(ylim[1], ylim[2], length.out=nc), 
                        ylim = ylim,
                        type="n", 
                        xaxs = "i", 
                        xaxt='n', yaxt='n', 
                        xlab="", ylab="")        
                                           
                   # user-defined calendar axis
                   if (is.null(date.format)) { time.tick = as.Date(spec.time.axis$at, tz=date.tz) }
                   if (!is.null(date.format)) { time.tick = as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz) }
                   
                   time.tick.label = spec.time.axis$labels
                   
                   # NAs among tick marks?
                   which.na = which(is.na(time.tick))
                   # NA warning
                   if ( length(which.na)>0 ) { time.axis.warning.na = T } 
                    
                   # in case of NAs among tick marks:
                   # NA tick marks (and corresponding labels) are omitted automatically by the axis command
                   
                   if ( is.logical(time.tick.label) ) {
                        # then time.tick.label has length 1 by criterion and IS TRUE OR FALSE
                        if (time.tick.label == T) {time.tick.label = time.tick}
                   }   
          
                   axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
                        las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                   
               } # end if (!time.axis.default)             
               
               mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)  
               
           }  # end if (show.date == T)  
           
       } # end if (label.time.axis == T)                            
  
    } # if (plot.rec == T)
    
    if (chronology.warning) { 
        warning("\nPlease check your calendar dates, format and time zone: dates may not be in an unambiguous format or chronological. The default numerical axis was used instead.") 
    }  
    
    if (time.axis.warning == T) {
        warning("\nPlease check your time axis specifications. Default settings were used.")
    }
    if (time.axis.warning.na == T) {
        warning("\nNAs were produced with your time axis specifications.")
    }
    
    ####################################
    ## Output
    ####################################
    
    output <- list(series = x.rec,
                 rec.waves = rec.waves,
                 loess.span = loess.span,
                 lvl = lvl,
                 only.coi = only.coi,
                 only.sig = only.sig, siglvl = siglvl, 
                 only.ridge = only.ridge,                
                 rnum.used = rnum.used,   
                 rescale = rescale,
                 dt = dt, dj = dj,
                 Period = Period, Scale = Scale,
                 nc = nc, nr = nr, 
                 axis.1 = WT$axis.1,
                 axis.2 = WT$axis.2,
                 date.format = date.format, date.tz = date.tz
                 )
                 
    class(output) = "reconstruct"

    out("Class attributes are accessible through following names:\n")
    out(names(output), "\n")             
     
    return(invisible(output))
       
      
}
