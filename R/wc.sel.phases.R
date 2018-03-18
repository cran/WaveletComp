wc.sel.phases <-
function(WC, sel.period = NULL, sel.lower = NULL, sel.upper = NULL, 
      only.coi = FALSE, 
      only.sig = TRUE, which.sig = "wp", siglvl = 0.05, 
      phase.cols = c("red", "blue"), 
      show.Angle = TRUE, use.sAngle = FALSE, Angle.col = "black", 
      show.legend = TRUE, legend.coords = "topleft", legend.horiz = TRUE,
      label.time.axis = TRUE, 
      show.date = FALSE, date.format = NULL, date.tz = NULL, 
      timelab = NULL, timetck = 0.02, timetcl = 0.5,
      spec.time.axis = list(at = NULL, labels = TRUE, 
                            las = 1, hadj = NA, padj = NA),
      label.phase.axis = TRUE, 
      phaselab = NULL, phasetck = 0.02, phasetcl = 0.5,
      spec.phase.axis = list(at = NULL, labels = TRUE, 
                             las = 1, hadj = NA, padj = NA),
      phaselim = c(-pi,pi+show.legend*ifelse(legend.horiz,0.8,2)),      
      main = NULL, sub = NULL,
      lwd = 1, lwd.Angle = 2, lwd.axis = 1, 
      verbose = FALSE) {      
                     
   ##########################################                  

   if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    } 
    
   ##########################################
  
   if (class(WC) == 'analyze.wavelet') { stop("Your object class is 'analyze.wavelet' --- please use wt.sel.phases!") }
  
   ##########################################
    
   nc = WC$nc    
   nr = WC$nr 
    
   series.data = WC$series
   
   # retrieve series names
   if (names(series.data)[1] == 'date') { my.series.names = names(series.data)[2:3] }
   if (names(series.data)[1] != 'date') { my.series.names = names(series.data)[1:2] }
   
   #######################################################################################
   ## Select periods
   #######################################################################################                  
                     
   Period = WC$Period

   # sel.period available?
   if (length(sel.period) != 0) {
   
       # nearest available period:
       sel.rnum = which(abs(Period-sel.period) == min(abs(Period-sel.period)))
       
   }    
   
   # in case, sel.period is not available, refer to sel.upper, sel.lower
   if (length(sel.period) == 0) {
   
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
       
   }    
       
   # selected band / range of periods   
   sel.period.band  = Period[sel.rnum]  
   sel.period.range = as.character(round(sel.period.band, 3))
   
   #######################################################################################
   ## Select corresponding phases and phase differences
   ####################################################################################### 
   
   Phase.x = WC$Phase.x
   Phase.y = WC$Phase.y
   Angle  = WC$Angle
   sAngle = WC$sAngle
   
   # only cone of influence? 
   if (only.coi == T) {
      
       for (i in (1:nc)) {
            for (s.ind in seq_len(nr)) {
                 if (WC$Scale[s.ind] > 2^WC$coi.2[i]) {
                     Phase.x[s.ind, i] = NA
                     Phase.y[s.ind, i] = NA
                     Angle[s.ind, i]   = NA
                     sAngle[s.ind, i]  = NA
                  }
             }
        } 
   }  
    
   # use significant parts only?
   if (only.sig == T) {
       
       if ( (which.sig == 'wp') & !is.null(WC$Power.xy.pval) ) { 
       
             ind.sig = (WC$Power.xy.pval < siglvl) 
             
             Phase.x[ind.sig == 0] = NA     
             Phase.y[ind.sig == 0] = NA
             Angle[ind.sig == 0] = NA
             sAngle[ind.sig == 0] = NA
            
       }
       if ( (which.sig == 'wc') & !is.null(WC$Coherence.pval) ) { 
       
             ind.sig = (WC$Coherence.pval < siglvl) 
             
             Phase.x[ind.sig == 0] = NA     
             Phase.y[ind.sig == 0] = NA
             Angle[ind.sig == 0] = NA
             sAngle[ind.sig == 0] = NA
             
       }
       if ( (which.sig == 'wt') & !is.null(WC$Power.x.pval) & !is.null(WC$Power.y.pval) ) { 
       
             ind.sig = (WC$Power.x.pval < siglvl) * (WC$Power.y.pval < siglvl) 
             ind.sig.x = (WC$Power.x.pval < siglvl)
             ind.sig.y = (WC$Power.y.pval < siglvl)
             
             Phase.x[ind.sig.x == 0] = NA     
             Phase.y[ind.sig.y == 0] = NA
             Angle[ind.sig == 0] = NA
             sAngle[ind.sig == 0] = NA
       }
            
   }  
   if (only.sig == F) { siglvl = NA } 
    
    
   # selected periods (period bands) 
   Phase.x = Phase.x[sel.rnum,]
   Phase.y = Phase.y[sel.rnum,]
   Angle  = Angle[sel.rnum,]
   sAngle = sAngle[sel.rnum,]
   
   # in case, several periods have been selected, compute column averages
   if (length(sel.rnum) > 1) {
       sel.period.range = paste(round(range(sel.period.band), 3), collapse=' - ')
      
           Phase.x = colMeans(Phase.x, na.rm=T)
           Phase.y = colMeans(Phase.y, na.rm=T)
           Angle   = colMeans(Angle, na.rm=T)
           sAngle  = colMeans(sAngle, na.rm=T)   
   }
   
   #######################################################################################
   ## Checks
   #######################################################################################
   
   # individual phase axis parameters
  
   if (!is.list(spec.phase.axis))       spec.phase.axis = list()
   if (is.null(spec.phase.axis$at))     spec.phase.axis$at = NULL
   if (is.null(spec.phase.axis$labels)) spec.phase.axis$labels = T
   if (is.null(spec.phase.axis$las))    spec.phase.axis$las = 1
   if (is.null(spec.phase.axis$hadj))   spec.phase.axis$hadj = NA
   if (is.null(spec.phase.axis$padj))   spec.phase.axis$padj = NA
  
   # initialize warning indicators for the case of individual phase axis specification
   # warning: reset to phase axis default?
   phase.axis.warning = F  
   # warning: NAs among phase axis tick marks?
   phase.axis.warning.na = F
  
   if ( (!is.null(spec.phase.axis$at)) & (label.phase.axis==F) )  warning("\nPlease set label.phase.axis = TRUE to make phase axis specification effective.", immediate. = TRUE)
   
   # date parameters
  
   if (is.null(date.format)) { date.format = WC$date.format }
   if (is.null(date.tz)) { date.tz = ifelse(is.null(WC$date.tz),"",WC$date.tz) }
  
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
   
   
   #######################################################################################
   ## start plotting
   #######################################################################################
      
   if (is.null(sub)) { sub=paste('selected period: ', sel.period.range, sep='') }  
   
#    matplot(WC$axis.1, data.frame(Phase.x, Phase.y), type = 'l', 
   matplot(1:nc, data.frame(Phase.x, Phase.y), type = 'l', 
           lty = 1, lwd = lwd, col = phase.cols, 
           xaxs = 'i', xaxt = 'n', 
           ylim = phaselim, yaxt = 'n',
           xlab = '', ylab = '', 
           main = main,
           sub = sub)
           
   legend.angle = NULL        
           
   if (show.Angle == T) {        
       if (use.sAngle == T) {        
#            lines(WC$axis.1, sAngle, lty=2, lwd=lwd.Angle, col=Angle.col)
           lines(1:nc, sAngle, lty=2, lwd=lwd.Angle, col=Angle.col)
           legend.angle = 'phase difference\n(smoothed)'
       } 
       if (use.sAngle == F) {        
#            lines(WC$axis.1, Angle, lty=2, lwd=lwd.Angle, col=Angle.col)
           lines(1:nc, Angle, lty=2, lwd=lwd.Angle, col=Angle.col)
           legend.angle = 'phase difference'
       } 
   }
   if (show.legend == T) {
       legend(legend.coords, legend=c(my.series.names, legend.angle), horiz=legend.horiz, col=c(phase.cols,ifelse(show.Angle,Angle.col,NULL)), 
              lty = c(1,1,ifelse(show.Angle,2,NULL)), lwd = c(lwd,lwd, ifelse(show.Angle,lwd.Angle,NULL)), 
              box.lwd=lwd.axis, 
              text.font = par()$font.lab, cex=par()$cex.lab)
   }
   
   ############################
   
   # label phase axis ?  
   if (label.phase.axis == T) {
   
       # check if there is a user-defined phase axis specification
       phase.axis.default = ( is.null(spec.phase.axis$at) )
      
       # first: check conditions of phase axis specification
       # phase.axis.warning = TRUE leads to default reset
      
       if ( !phase.axis.default ) {
      
            # 1. Are tick marks (spec.phase.axis$at) interpretable as numbers?
            # is.numeric returns TRUE iff there is no character value but at least one non-NA value which is interpretable as number
            # (further logical entries T and F are interpreted as 1 and 0 which could be accepted here...)
           
            # if numeric: are phase values within appropriate intervall and all finite? (NAs, being also interpreted as not finite, will be removed during plotting)
            if ( is.numeric(spec.phase.axis$at) ) { 
                 phase.axis.warning = ( sum(!is.na(spec.phase.axis$at)) != sum(is.finite(spec.phase.axis$at)) ) 
            }
            else { phase.axis.warning = TRUE }
           
            # 2. Are labels (spec.phase.axis$labels) appropriate? Either logical and single-valued or non-logical and of the same length as tick marks?
            # is.logical returns TRUE iff there is not any character or numeric value, but there could be NA, even if it is the only value.
           
            # if logical:  single-valued and non-NA? 
            if ( is.logical(spec.phase.axis$labels) ) { 
                 phase.axis.warning = ( phase.axis.warning | ifelse(length(spec.phase.axis$labels)!=1,TRUE, is.na(spec.phase.axis$labels)) )
            }
           
            # if non-logical: do vectors of tick marks and labels have equal length?
            if ( !is.logical(spec.phase.axis$labels) ) {
                phase.axis.warning = ( phase.axis.warning | (length(spec.phase.axis$labels) != length(spec.phase.axis$at)) ) 
            }
           
       } # end if ( !phase.axis.default )
      
       # default reset in case of any warning
       phase.axis.default = (phase.axis.default | phase.axis.warning)
   
       if ( (is.null(phaselab)) | (!is.null(phaselab) & phase.axis.warning) ) {
           phaselab='phase'
           if (length(sel.rnum) > 1) {
               phaselab='average phase'
           }
        }
       
        if (phase.axis.default) {
        
            axis(2, lwd = lwd.axis, at=seq(-pi, pi, pi/2), labels=NA, tck=phasetck, tcl=phasetcl)
            axis(4, lwd = lwd.axis, at=seq(-pi, pi, pi/2), labels=NA, tck=phasetck, tcl=phasetcl)
            mtext(c(expression(-pi),expression(-pi/2),0,expression(pi/2),expression(pi)), at=seq(-pi, pi, pi/2), side = 2, line = par()$mgp[2]-0.5, las=1, 
                  font = par()$font.axis, cex=par()$cex.axis)
        }
        
        if (!phase.axis.default) {
      
            phase.tick = spec.phase.axis$at
            phase.tick.label = spec.phase.axis$labels
          
            # NAs among tick marks?
            which.na = which(is.na(phase.tick))
            # NA warning
            if ( length(which.na)>0 ) { phase.axis.warning.na = T } 
           
            # in case of NAs among tick marks:
            # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
            axis(2, lwd = lwd.axis, at = phase.tick, labels = phase.tick.label, tck=phasetck, tcl=phasetcl,
                 las = spec.phase.axis$las, hadj = spec.phase.axis$hadj, padj=spec.phase.axis$padj,  
                 mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
            axis(4, lwd = lwd.axis, at = phase.tick, labels = NA, tck=phasetck, tcl=phasetcl)    

        } # end if (!phase.axis.default)
        
        mtext(phaselab, side = 2, line = par()$mgp[1]-0.5, font = par()$font.lab, cex=par()$cex.lab)
       
   } # end  if (label.phase.axis == T)
   
   ############################
   
   my.date = NULL
   
   #########################
   
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
        
      ##  
  
      if ( (is.null(timelab) & time.axis.default) | (!is.null(timelab) & time.axis.warning) )  {timelab=ifelse(show.date, 'calendar date','index')}
        
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
                   plot(my.date, rep(0,nc), 
                        ylim=phaselim,
                        type="n", 
                        xaxs = "i", 
                        yaxt='n', 
                        xlab="", ylab="",
                        lwd = lwd.axis, tck=timetck, tcl=timetcl, 
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                        
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   # empty plot, but calendar axis (w/o plotting)
                   plot(my.date, rep(0,nc), 
                        ylim=phaselim, 
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
      
   }  # end if (label.time.axis == T)     
  
  ####################################
  ## Output
  ####################################
  
  if (phase.axis.warning == T) {
      warning("\nPlease check your phase axis specifications. Default settings were used.")
  }
  if (phase.axis.warning.na == T) {
      warning("\nNAs were produced with your phase axis specifications.")
  }
  
  if (chronology.warning) { 
      warning("\nPlease check your calendar dates, format and time zone: dates may not be in an unambiguous format or chronological. The default numerical axis was used instead.") 
  }  
  
  if (time.axis.warning == T) {
      warning("\nPlease check your time axis specifications. Default settings were used.")
  }
  if (time.axis.warning.na == T) {
      warning("\nNAs were produced with your time axis specifications.")
  }
    
  output <- list(Period = sel.period.band,
                 Phase.x = Phase.x,
                 Phase.y = Phase.y,
                 Angle = Angle,
                 sAngle = sAngle,
                 only.coi = only.coi,
                 only.sig = only.sig, which.sig = which.sig, siglvl = siglvl, 
                 date = my.date, date.format = date.format, date.tz = date.tz,
                 axis.1 = WC$axis.1
                 )
                 
  class(output) = "sel.phases"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")             
     
  return(invisible(output))
         

}
