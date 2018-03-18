
### DEFINITION:
### (v2018-02-14)

periodic.series = function(start.period = 100, end.period = start.period, phase = 0, length = 600, make.plot = FALSE)
{
  if (start.period < end.period) { m = start.period/end.period }
  if (start.period >= end.period) { m = end.period/start.period }
  N = length/(1-0.5*(1-m))
  x = 0:(N-1)
  if (start.period < end.period) { period = end.period }
  if (start.period >= end.period) { period = start.period }
  p = rep(period, length = N)
  y = sin((x + phase) * 2 * pi/p)
  if (start.period < end.period) { y = y[N:1] }
  # "thinning out":
  i = 1:N
  cum.keep.share = 1 - 0.5*i*(1-m)/N
  keep = floor(cum.keep.share*i)
  diff.keep = c(1, diff(keep))
  if (sum(diff.keep) < length) {diff.keep[match(0, diff.keep)] = 1}
  if (sum(diff.keep) > length) {diff.keep[match(1, diff.keep)] = 0}
  if (start.period < end.period) {
    # modify diff.keep:
    first.1.rev = match(1,diff.keep[N:1]) # first one in reversed diff.keep
    if (first.1.rev > 1)  {
      # cut off 0s at the end, add at beginning:
      diff.keep = c(rep(0, first.1.rev - 1), diff.keep[1:(N-(first.1.rev - 1))])
    }
  }
  if (start.period >= end.period) {
    # modify diff.keep:
    first.1 = match(1,diff.keep) # first one in diff.keep
    if (first.1 > 1)  {
      # cut off 0s at the beginning, add at the end:
      diff.keep = c(diff.keep[first.1:N], rep(0, first.1 - 1))
    }
  }
  z = y[diff.keep == 1]
  if (start.period < end.period) {z = rev(z)}
  if (make.plot) {ts.plot(z, xlab = 'time', ylab = 'periodic series')}
  return(z)
}


### EXAMPLE:

# x1 = 0.8*periodic.series(start.period = 100, end.period = 95, phase = 0, length = 1000)
# x2 =     periodic.series(start.period = 100, end.period = 100, phase = 0, length = 1000)
# x3 = 1.2*periodic.series(start.period = 100, end.period = 105, phase = 0, length = 1000)
# 
# ts.plot(x2, ylim = c(-2, +2), xlab = 'time', ylab = 'series with variable period')
# lines(x1, col = 'blue')
# lines(x3, col = 'red')
# legend('topleft', legend = c('speeding up (end period = 95)', 'period = 100', 'slowing down (end period = 105)'), lty = 1, col = c('blue', 'black', 'red'))
