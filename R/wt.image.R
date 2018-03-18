wt.image <-
function(WT, my.series = 1, exponent = 1,
     plot.coi = TRUE, 
     plot.contour = TRUE, siglvl = 0.1, col.contour = "white",
     plot.ridge = TRUE, lvl = 0, col.ridge = "black", 
     color.key = "quantile", 
     n.levels = 100, color.palette = "rainbow(n.levels, start = 0, end = .7)",
     maximum.level = NULL,
     useRaster = TRUE, max.contour.segments = 250000,
     plot.legend = TRUE,
     legend.params = list(width = 1.2, shrink = 0.9, mar = 5.1, 
                          n.ticks = 6, 
                          label.digits = 1, label.format = "f", 
                          lab = NULL, lab.line = 2.5),
     label.time.axis = TRUE, 
     show.date = FALSE, date.format = NULL, date.tz = NULL, 
     timelab = NULL, timetck = 0.02, timetcl = 0.5,
     spec.time.axis = list(at = NULL, labels = TRUE, 
                           las = 1, hadj = NA, padj = NA),
     label.period.axis = TRUE, 
     periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
     spec.period.axis = list(at = NULL, labels = TRUE, 
                             las = 1, hadj = NA, padj = NA),
     main = NULL,
     lwd = 2, lwd.axis = 1, 
     graphics.reset = TRUE,
     verbose = FALSE) {       
                    
  ################################################    
                     
  if(verbose == T){
       out <- function(...){ cat(...) }
    }
    else{
       out <- function(...) { }
    }  

  default.options = options() 
  
  options(max.contour.segments = as.integer(max.contour.segments))
   
  #################################################
  
  axis.1 <- WT$axis.1
  axis.2 <- WT$axis.2 
   
  series.data = WT$series
  
  if (exponent <= 0) { stop("Please use a positive exponent, or return to default setting (exponent=1)!") }  
    
  ####################################
  ## Identify the scenario
  ####################################
    
  if (class(WT) == 'analyze.wavelet') {
    
      out("Your input object class is 'analyze.wavelet'...\n") 
         
      my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[2], names(series.data)[1]) 
        
      Power = WT$Power^exponent
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
          Power = WT$Power.x^exponent
          Power.pval = WT$Power.x.pval
          Ridge = WT$Ridge.x
      }
      if (which.series.num == 2) {
          Power = WT$Power.y^exponent
          Power.pval = WT$Power.y.pval
          Ridge = WT$Ridge.y
      }      
        
  }   
  
  out(paste("A wavelet power image of your time series '", my.series, "' will be plotted...", sep=''),'\n') 
  
  if ( (!is.null(maximum.level)) & (is.element(color.key, c('quantile','q'))) ) { 
       warning("\nPlease set color.key = 'i' to make your maximum level specification effective.", immediate. = TRUE) 
  }
  if ( is.element(color.key, c('interval','i')) ) {  
       if (is.null(maximum.level)) {maximum.level = max(Power)} 
       if (maximum.level < max(Power)) { 
           stop(paste("... plot can't be produced! Your choice of maximum plot level is smaller than the maximum level observed! Please choose maximum.level larger than ", max(Power), " or return to default setting (maximum.level = NULL)!", sep='')) 
       }
  }
           
  if (is.element(color.key,c('interval','i'))) {    
      wavelet.levels = seq(from=0, to=maximum.level, length.out=n.levels+1)
  }  
  if (is.element(color.key,c('quantile','q'))) {  
      wavelet.levels = quantile(Power, probs=seq(from=0, to=1, length.out=n.levels+1)) 
  }   
  key.cols = rev(eval(parse(text=color.palette)))
  
  # legend parameters  
  
  if (!is.list(legend.params))              legend.params = list()
  if (is.null(legend.params$width))         legend.params$width = 1.2
  if (is.null(legend.params$shrink))        legend.params$shrink = 0.9
  if (is.null(legend.params$mar)) legend.params$mar = ifelse(is.null(legend.params$lab), 5.1, 6.1)   
  if (is.null(legend.params$n.ticks))       legend.params$n.ticks = 6
  if (is.null(legend.params$label.digits))  legend.params$label.digits = 1
  if (is.null(legend.params$label.format))  legend.params$label.format = "f"
  if (is.null(legend.params$lab.line))      legend.params$lab.line = 2.5
  
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
  
  # individual period axis parameters
  
  if (!is.list(spec.period.axis))       spec.period.axis = list()
  if (is.null(spec.period.axis$at))     spec.period.axis$at = NULL
  if (is.null(spec.period.axis$labels)) spec.period.axis$labels = T
  if (is.null(spec.period.axis$las))    spec.period.axis$las = 1
  if (is.null(spec.period.axis$hadj))   spec.period.axis$hadj = NA
  if (is.null(spec.period.axis$padj))   spec.period.axis$padj = NA
  
  # initialize warning indicators for the case of individual period axis specification
  # warning: reset to period axis default?
  period.axis.warning = F  
  # warning: NAs among period axis tick marks?
  period.axis.warning.na = F
    
  if ( (!is.null(spec.period.axis$at)) & (label.period.axis==F) )  warning("\nPlease set label.period.axis = TRUE to make period axis specification effective.", immediate. = TRUE)
  
  #######################################################################################
  ## start plotting
  #######################################################################################
  
  op = par(no.readonly = TRUE)
  
  image.plt  = par()$plt
  legend.plt = NULL
  
  if (plot.legend == T) {
  
      # construct plot regions for image and legend
   
      legend.plt = par()$plt
      
      char.size = par()$cin[1]/par()$din[1]
    
      hoffset       = char.size * par()$mar[4]
      legend.width  = char.size * legend.params$width
      legend.mar    = char.size * legend.params$mar
      
      legend.plt[2] = 1 - legend.mar
      legend.plt[1] = legend.plt[2] - legend.width
  
      vmar = (legend.plt[4] - legend.plt[3]) * ((1 - legend.params$shrink)/2)
  
      legend.plt[4] = legend.plt[4] - vmar
      legend.plt[3] = legend.plt[3] + vmar

      image.plt[2] = min(image.plt[2], legend.plt[1] - hoffset)
      
      # plot legend first
               
      par(plt = legend.plt)

      key.marks  = round(seq(from = 0, to = 1, length.out=legend.params$n.ticks)*n.levels)
      key.labels = formatC(as.numeric(wavelet.levels), digits = legend.params$label.digits, format = legend.params$label.format)[key.marks+1]

      image(1, seq(from = 0, to = n.levels), matrix(wavelet.levels, nrow=1), col = key.cols, breaks = wavelet.levels, useRaster=T, xaxt='n', yaxt='n', xlab='', ylab='')
      axis(4, lwd=lwd.axis, at=key.marks, labels=NA, tck=0.02, tcl=(par()$usr[2]-par()$usr[1])*legend.params$width-0.04)
      mtext(key.labels, side = 4, at = key.marks, line = 0.5, las=2, font = par()$font.axis, cex=par()$cex.axis)
      text(x = par()$usr[2] + (1.5+legend.params$lab.line)*par()$cxy[1], y=n.levels/2, labels=legend.params$lab, xpd=NA, srt = 270, font = par()$font.lab, cex=par()$cex.lab)
  
      box(lwd = lwd.axis)    
    
      par(new=TRUE, plt = image.plt)  
      
  }    
  
  #######################################################################################
  ## plot power image
  #######################################################################################
  
  image(axis.1, axis.2, t(Power), col = key.cols, breaks = wavelet.levels, useRaster = useRaster,
        ylab = "", xlab = '', axes = FALSE, main = main)        
     
  # plot contour lines?     
  if  ((plot.contour == T) & (is.null(Power.pval) == F)) {      
      contour(axis.1, axis.2, t(Power.pval) < siglvl, levels = 1, lwd = lwd, 
              add = TRUE, col = col.contour, drawlabels = FALSE)
  }
  
  # plot ridge?
  if  (plot.ridge == T) {    
       Ridge = Ridge * (Power >= lvl)      
       contour(axis.1, axis.2, t(Ridge), levels = 1, lwd = lwd,
              add = TRUE, col = col.ridge, drawlabels = FALSE)         
  }     

  # plot cone of influence?
  if (plot.coi == T) {
      polygon(WT$coi.1, WT$coi.2, border = NA, col = rgb(1, 1, 1, 0.5))
  }  
  
  box(lwd = lwd.axis)
  
  ##############################################
  
  # label period axis ?  
  if (label.period.axis == T) {
  
      # check if there is a user-defined period axis specification
      period.axis.default = ( is.null(spec.period.axis$at) )
      
      # first: check conditions of period axis specification
      # period.axis.warning = TRUE leads to default reset
      
      if ( !period.axis.default ) {
      
           # 1. Are tick marks (spec.period.axis$at) interpretable as numbers?
           # is.numeric returns TRUE iff there is no character value but at least one non-NA value which is interpretable as number
           # (further logical entries T and F are interpreted as 1 and 0 which could be accepted here... entry 0 will be removed later on)
           
           # if numeric: are period values all positive? all finite? (NAs, being also interpreted as not finite, will be removed during plotting)
           if ( is.numeric(spec.period.axis$at) ) { 
                period.axis.warning = ( ( sum(spec.period.axis$at<=0, na.rm=T)>0 ) | ( sum(!is.na(spec.period.axis$at)) != sum(is.finite(spec.period.axis$at)) ) )
           }
           else { period.axis.warning = TRUE }
           
           # 2. Are labels (spec.period.axis$labels) appropriate? Either logical and single-valued or non-logical and of the same length as tick marks?
           # is.logical returns TRUE iff there is not any character or numeric value, but there could be NA, even if it is the only value.
           
           # if logical:  single-valued and non-NA? 
           if ( is.logical(spec.period.axis$labels) ) { 
                period.axis.warning = ( period.axis.warning | ifelse(length(spec.period.axis$labels)!=1,TRUE, is.na(spec.period.axis$labels)) )
           }
           
           # if non-logical: do vectors of tick marks and labels have equal length?
           if ( !is.logical(spec.period.axis$labels) ) {
               period.axis.warning = ( period.axis.warning | (length(spec.period.axis$labels) != length(spec.period.axis$at)) ) 
           }
           
      } # end if ( !period.axis.default )
      
      # default reset in case of any warning
      period.axis.default = (period.axis.default | period.axis.warning)
  
      if ( (is.null(periodlab)) | (!is.null(periodlab) & period.axis.warning) ) {periodlab='period'}
    
      if (period.axis.default) {
          period.tick = unique(trunc(axis.2))
          period.tick[period.tick<log2(WT$Period[1])] = NA
          period.tick = na.omit(period.tick)
          period.tick.label = 2^(period.tick)   
          axis(2, lwd = lwd.axis, at = period.tick, labels = NA, tck=periodtck, tcl=periodtcl)
          axis(4, lwd = lwd.axis, at = period.tick, labels = NA, tck=periodtck, tcl=periodtcl)
          mtext(period.tick.label, side = 2, at = period.tick, las = 1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
      }
      
      if (!period.axis.default) {
      
          period.tick = log2(spec.period.axis$at)
          period.tick[(period.tick<log2(WT$Period[1]))] = NA
          # (works even if any period.tick value is NA beforehand)   
          period.tick.label = spec.period.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(period.tick))
          # NA warning
          if ( length(which.na)>0 ) { period.axis.warning.na = T } 
           
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          if ( is.logical(period.tick.label) ) {
               # then period.tick.label has length 1 by criterion
               # The following works even if any period.tick value is NA   
               if (period.tick.label == T) { period.tick.label = 2^(period.tick) }
          }     
          
          axis(2, lwd = lwd.axis, at = period.tick, labels = period.tick.label, tck=periodtck, tcl=periodtcl,
               las = spec.period.axis$las, hadj = spec.period.axis$hadj, padj=spec.period.axis$padj, 
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          axis(4, lwd = lwd.axis, at = period.tick, labels = NA, tck=periodtck, tcl=periodtcl)    

      } # end if (!period.axis.default)
      
      mtext(periodlab, side = 2, line = par()$mgp[1]-0.5, font = par()$font.lab, cex=par()$cex.lab)
          
  }
  
  #################################
   
  # label time axis ? --- checks whether time.axis.default or not is applicable (including format and chronology tests of date in case show.date=T
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
      
  }  # end if (label.time.axis == T)   
   
  ##############################################
  
  half.increment.2 = (axis.2[2]-axis.2[1])/2
  my.ylim = c(min(axis.2)-half.increment.2, max(axis.2)+half.increment.2)
  
  index.condition = ( (show.date==F) | (label.time.axis==F) )
  
  if ((WT$dt != 1) & (index.condition)) {
  
  par(new=TRUE)
  
  # empty plot, but index axis (w/o plotting)
  plot(1:WT$nc, seq(min(axis.2),max(axis.2), length.out=WT$nc),
       xlim = c(0.5,WT$nc+0.5), 
       ylim = my.ylim, 
       type="n", 
       xaxs = "i", yaxs ='i', 
       xaxt='n', yaxt='n', 
       xlab="", ylab="")
       
  }     
       
  ##############################################
  
  # finally: label time axis?  
  if (label.time.axis == T) { 
    
      if ( (is.null(timelab) & time.axis.default) | (!is.null(timelab) & time.axis.warning) ) {timelab=ifelse(show.date, 'calendar date','index')}
      
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
      
               half.increment = difftime(my.date[2],my.date[1], tz=date.tz)/2
                     
               if (is.null(date.format)) { my.xlim = as.Date(range(my.date) - c(half.increment,-half.increment),tz=date.tz) } 
               else { my.xlim = as.POSIXct(range(my.date) - c(half.increment,-half.increment),format=date.format, tz=date.tz) }
                 
               par(new=TRUE)
               
               if (time.axis.default) {
               
                   # empty plot, but calendar axis (w/ plotting)
                   plot(my.date, seq(min(axis.2),max(axis.2), length.out=WT$nc), 
                        xlim = my.xlim, 
                        ylim = my.ylim, 
                        type="n", 
                        xaxs = "i", yaxs ='i', 
                        yaxt='n', 
                        xlab="", ylab="",
                        lwd = lwd.axis, tck=timetck, tcl=timetcl, 
                        mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
                        
               } # end if (time.axis.default)
               
               if (!time.axis.default) {
               
                   # empty plot, but calendar axis (w/o plotting)
                   plot(my.date, seq(min(axis.2),max(axis.2), length.out=WT$nc),
                        xlim = my.xlim, 
                        ylim = my.ylim, 
                        type="n", 
                        xaxs = "i", yaxs ='i', 
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
      
  }  # end finally: if (label.time.axis == T)

  
  #######################################################################################
  ## apropos graphical parameters
  #######################################################################################
  
  # reset contour line options
  options(default.options)
  
  # reset graphical parameters?
  if (graphics.reset == T) {
      par(op)
  }     
  
  if (period.axis.warning == T) {
      warning("\nPlease check your period axis specifications. Default settings were used.")
  }
  if (period.axis.warning.na == T) {
      warning("\nNAs were produced with your period axis specifications.")
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
  
  # output of graphical parameters
  
  output = list(op = op, image.plt = image.plt, legend.plt=legend.plt)
  class(output) = "graphical parameters"

  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")    
    
  return(invisible(output)) 
}
