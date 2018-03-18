## Version 1.1 (18 March 2018)
Tools for displaying and analyzing periodic phenomena across time have been extended. The main innovations are:
  * All functions of family wt.<> (showing results concerning a single time series) can now also be applied to extract univariate outcomes from cross-wavelet and coherence analysis (objects of class "analyze.coherency").
  * It is possible to control the color gradation of time-period spectrum plots, and accentuate the contrast, by raising the wavelet power values to any (positive) exponent before plotting.
  * Setting a maximum level for the color bar facilitates the visual comparison of time-period spectrum plots. Maximum and minimum plot levels are options for plots of averages too.
  * The time and period axes are now easier to individualize by specifying tick marks and labels. Coordinates on the time axis can be conveniently addressed via an index or a POSIXct object.
  * Graphical parameters of global coverage (cex.axis, font.axis, cex.lab, font.lab, mgp etc., see par) as well as parameters of local coverage (within axis specification options) help fine-tune plots.
  * Two more real-world data sets have been included in WaveletComp, namely:
    o Data set weather.radiation.Mannheim, containing daily weather and ambient radiation readings from Mannheim (Germany).
    o Data set USelection2016.Instagram, containing hourly numbers of candidate-related media uploads to Instagram right before the 2016 US presidential election.

