wt.avg <-
function(WT, my.series = 1, exponent = 1,
  show.siglvl = TRUE, 
  siglvl = c(0.05, 0.1), sigcol = c("red", "blue"), sigpch = 20, sigcex = 1,
  minimum.level = NULL, maximum.level = NULL,
  label.avg.axis = TRUE, 
  averagelab = NULL, averagetck = 0.02, averagetcl = 0.5,
  spec.avg.axis = list(at = NULL, labels = TRUE, 
                       las = 1, hadj = NA, padj = NA),
  label.period.axis = TRUE, 
  periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
  spec.period.axis = list(at = NULL, labels = TRUE, 
                          las = 1, hadj = NA, padj = NA),
  show.legend = TRUE, legend.coords = "topright", 
  main = NULL, 
  lwd = 1, col = 1, 
  lwd.axis = 1,
  verbose = FALSE) {       
         
  ##########################################       
 
  if(verbose == T){
       out <- function(...){ cat(...) }
    }
  else{
       out <- function(...) { }
    }  
           
  if (exponent <= 0) { stop("Please use a positive exponent, or return to default setting (exponent=1)!") }         
         
  ####################################
  ## Identify the scenario
  ####################################
  
  series.data = WT$series
    
  if (class(WT) == 'analyze.wavelet') {
    
      out("Your input object class is 'analyze.wavelet'...\n")       
      
      my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[2], names(series.data)[1])
              
      Power.avg = WT$Power.avg^exponent
      Power.avg.pval = WT$Power.avg.pval   
        
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
          
          Power.avg = WT$Power.x.avg^exponent
          Power.avg.pval = WT$Power.x.avg.pval
          
      }
      if (which.series.num == 2) {
         
          Power.avg = WT$Power.y.avg^exponent
          Power.avg.pval = WT$Power.y.avg.pval
          
      }      
        
  }  
  
  out(paste("Power averages across time of your time series '", my.series, "' will be plotted...", sep=''), '\n') 
  
  if (is.null(minimum.level)) {minimum.level = min(Power.avg)}
  if (minimum.level > min(Power.avg)) { 
      stop(paste("... no plot produced! Your choice of minimum plot level is larger than the minimum level observed! Please choose minimum.level smaller than ", min(Power.avg), " or return to default setting (minimum.level = NULL)!", sep='')) 
  } 
  if (is.null(maximum.level)) {maximum.level = max(Power.avg)}
  if (maximum.level < max(Power.avg)) { 
      stop(paste("... no plot produced! Your choice of maximum plot level is smaller than the maximum level observed! Please choose maximum.level larger than ", max(Power.avg), " or return to default setting (maximum.level = NULL)!", sep='')) 
  } 
  
  
  # individual average axis parameters
  
  if (!is.list(spec.avg.axis))       spec.avg.axis = list()
  if (is.null(spec.avg.axis$at))     spec.avg.axis$at = NULL
  if (is.null(spec.avg.axis$labels)) spec.avg.axis$labels = T
  if (is.null(spec.avg.axis$las))    spec.avg.axis$las = 1
  if (is.null(spec.avg.axis$hadj))   spec.avg.axis$hadj = NA
  if (is.null(spec.avg.axis$padj))   spec.avg.axis$padj = NA
  
  # initialize warning indicators for the case of individual average axis specification
  # warning: reset to average axis default?
  avg.axis.warning = F  
  # warning: NAs among average axis tick marks?
  avg.axis.warning.na = F
  
  if ( (!is.null(spec.avg.axis$at)) & (label.avg.axis==F) )  warning("\nPlease set label.avg.axis = TRUE to make average axis specification effective.", immediate. = TRUE)   
  
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

  plot(Power.avg, log2(WT$Period), xlim = range(c(Power.avg,minimum.level,maximum.level)), 
       lwd = lwd, col = col, type = "l",  
       axes = FALSE, 
       ylab = "", xlab = '', 
       yaxs = 'i', 
       main = main)
  
  
  if ((show.siglvl == T) & (is.null(Power.avg.pval) == F)) { 
     
     P.dat = data.frame(Pvalue = Power.avg.pval, Log.period = log2(WT$Period), Average = Power.avg)
     
     n.sig = length(siglvl)
     if (length(sigpch) == 1) {sigpch = rep(sigpch,n.sig)}
     if (length(sigcex) == 1) {sigcex = rep(sigcex,n.sig)}
     if (n.sig != length(sigpch)) {sigpch = rep(20,n.sig)}
     if (n.sig != length(sigcol)) {sigcol = 1:n.sig}
     if (n.sig != length(sigcex)) {sigcex = rep(1,n.sig)}
     
     siglvl.order = order(siglvl,decreasing=T)
     sig.params = data.frame(siglvl = siglvl[siglvl.order], sigcol=sigcol[siglvl.order], sigpch = sigpch[siglvl.order], sigcex = sigcex[siglvl.order])
     for (i in (1:length(siglvl))) {
         with(P.dat[P.dat$Pvalue < sig.params$siglvl[i],], points(Average, Log.period, pch = sig.params$sigpch[i], col = as.character(sig.params$sigcol[i]), cex = sig.params$sigcex[i]))
     } 
     
     
     if (show.legend == T) {
         legend(legend.coords, legend=siglvl, pch=sigpch, col=sigcol, horiz=F, box.lwd=lwd.axis, text.font = par()$font.lab, cex=par()$cex.lab, pt.cex=sigcex)
     }
      
  }
  
  box(lwd = lwd.axis)
  
  # label average axis ?  
  if (label.avg.axis == T) {
  
    # check if there is a user-defined average axis specification
      avg.axis.default = ( is.null(spec.avg.axis$at) )
      
      # first: check conditions of average axis specification
      # avg.axis.warning = TRUE leads to default reset
      
      if ( !avg.axis.default ) {
      
           # 1. Are tick marks (spec.avg.axis$at) interpretable as numbers?
           # is.numeric returns TRUE iff there is no character value but at least one non-NA value which is interpretable as number
           # (further logical entries T and F are interpreted as 1 and 0 which could be accepted here... entry 0 will be removed later on)
           
           # if numeric: are average values all non-negative? all finite? (NAs, being also interpreted as not finite, will be removed during plotting)
           if ( is.numeric(spec.avg.axis$at) ) { 
                avg.axis.warning = ( sum(!is.na(spec.avg.axis$at)) != sum(is.finite(spec.avg.axis$at)) ) 
           }
           else { avg.axis.warning = TRUE }
           
           # 2. Are labels (spec.avg.axis$labels) appropriate? Either logical and single-valued or non-logical and of the same length as tick marks?
           # is.logical returns TRUE iff there is not any character or numeric value, but there could be NA, even if it is the only value.
           
           # if logical:  single-valued and non-NA? 
           if ( is.logical(spec.avg.axis$labels) ) { 
                avg.axis.warning = ( avg.axis.warning | ifelse(length(spec.avg.axis$labels)!=1,TRUE, is.na(spec.avg.axis$labels)) )
           }
           
           # if non-logical: do vectors of tick marks and labels have equal length?
           if ( !is.logical(spec.avg.axis$labels) ) {
               avg.axis.warning = ( avg.axis.warning | (length(spec.avg.axis$labels) != length(spec.avg.axis$at)) ) 
           }
           
      } # end if ( !avg.axis.default )
      
      # default reset in case of any warning
      avg.axis.default = (avg.axis.default | avg.axis.warning)
  
      if ( (is.null(averagelab)) | (!is.null(averagelab) & avg.axis.warning) ) {averagelab='average wavelet power'}
   
      if (avg.axis.default) {
          A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = averagetck, tcl = averagetcl)
          mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
      }
      
      if (!avg.axis.default) {
      
          avg.tick = spec.avg.axis$at
          avg.tick.label = spec.avg.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(avg.tick))
          # NA warning
          if ( length(which.na)>0 ) { avg.axis.warning.na = T } 
                  
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          axis(1, lwd = lwd.axis, at = avg.tick, labels = avg.tick.label, tck = averagetck, tcl = averagetcl,
                las = spec.avg.axis$las, hadj = spec.avg.axis$hadj, padj=spec.avg.axis$padj,  
                mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)        

      } # end if (!avg.axis.default)
      
      mtext(averagelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)  
     
  } # end if (label.avg.axis == T)

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
          period.tick <- unique(trunc(WT$axis.2))
          period.tick[period.tick<log2(WT$Period[1])] = NA
          period.tick = na.omit(period.tick)
          period.tick.label <- 2^(period.tick)   
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
  
  } # end if (label.period.axis == T)
  
  
  if (avg.axis.warning == T) {
      warning("\nPlease check your specifications for the axis of averages. Default settings were used.")
  }
  if (avg.axis.warning.na == T) {
      warning("\nNAs were produced with your specifications for the axis of averages.")
  }
  
  if (period.axis.warning == T) {
      warning("\nPlease check your period axis specifications. Default settings were used.")
  }
  if (period.axis.warning.na == T) {
      warning("\nNAs were produced with your period axis specifications.")
  }
 
}
