rm( list=ls() )

### Calculates the fft of the fluctuations and fit Heilrich model (similar to Almonacid et al. 2018)

## folder to process in which the data are
fold = "MDA231"

library( xlsx )
library( ggplot2 )
library( stats )

## free memory
library( rJava )
jgc <- function()
{
  .jcall("java/lang/System", method = "gc")
}    

fluct = function( row )
{
    # Calculate the fluctuations around the mean radius
    return( (row - mean( row )) )
}

getFFT = function( data )
{
  # Calculate the fft of the fluctuations
  fftres = ( fft( data ) )
  ## normalise the fft by number of modes 
  fftres = fftres / (length(fftres[1,]))
  ## take the modulus of the fft
  fftres = apply( fftres, 2, abs )
  ## sqaure
  fftres = fftres^2
  ## keep only half of the mode (symetric)
  return( fftres[,1:(length(fftres[1,])/2)] )
}

getFFT2 = function( data )
{
  # Calculate the fft of the fluctuations
  fftres = ( fft( data ) )
  ## normalise the fft by number of modes 
  fftres = fftres / (length(fftres))
  ## take the modulus of the fft
  fftres = abs( fftres )
  ## square of the modulus
  fftres = fftres^2
  ## keep only half of the mode (symetric)
  return( fftres[1:(length(fftres)/2)] )
}

# Load the data
file_list = list.files( fold, pattern = "*.xlsx", full.names = TRUE )
options(java.parameters = "-Xmx8000m")

ffts = c()
for (file in file_list)
{
  gc()
  jgc()
  data = read.xlsx( file, sheetIndex=1, header=T )
  moy = data$Moyenne
  thetas = data$thetha
  data$Number = NULL
  data$thetha = NULL
  data$Moyenne = NULL
  fluctuations = apply( data, 1, fluct )
  ## check mean fluctuations squared
  print(mean(fluctuations^2))
  for (t in 1:length(fluctuations[,1]))
  {
    fft_cur = getFFT2( fluctuations[t,] )
    if (t==1)
    {
      fft_mat = fft_cur
    } else {
      fft_mat = cbind( fft_mat, fft_cur )
    }
  }
  ## average all the ffts in time
  fft_mean = apply( fft_mat, 1, mean )
  #print(length(fft_mean))
  #plot( (fft_mean)^2, type="l" )  
  ffts = cbind( ffts, fft_mean )
}

## plot to check visually
plot( ffts[,1], type="l" )
for (ind in 2:length(ffts[1,]))
{
  points( ffts[,ind], type="l" )
}

## write table results
ffts = ffts[2:length(ffts[,1]),]  ## remove the mode 0
ffts = cbind( seq(1:length(ffts[,1])), ffts )
colnames(ffts) = c( "Mode", file_list )
write.table( ffts, paste(fold, "results_fft.csv", sep="/"), sep=";", row.names=FALSE, col.names=TRUE )
