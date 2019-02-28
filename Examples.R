# Alpha-shape Fitting to Spectrum algorithm (AFS) &
# Alpha-shape and Lab Source Fitting to Spectrum algorithm (ALSFS)

source("functions/AFS.R")
source("functions/ALSFS.R")
source("functions/LS_Smoothing.R")

######################################################
# An example spectrum to show the effect of AFS and ALSFS.
# Read in the spectrum into variable x.
x <- read.csv("examples/ExampleSpectrum.csv")
plot(x, type = "l")

# Use q = 0.95 and d = 0.25 to run AFS function on the spectrum.
res1 <- AFS(x, .95, .25)
plot(x$wv, res1, type = "l")
abline(h = 1, col = "red")

# Read in the lab source spectrum to variable ls.
# Use q = 0.95 and d = 0.25 to run ALSFS function.
ls <- read.csv("examples/LabSource.csv")
plot(ls, type="l")
res2 <- ALSFS(x, ls, .95, .25)
plot(x$wv, res2, type = "l")
abline(h = 1, col = "red")


######################################################
# An exmaple to use AFS to smooth the raw lab source spectrum.
# Read in the raw lab source spectrum to variable rls.
rls <- read.csv("examples/RawLabSource.csv")
plot(rls, type="l")
res3 <- LSS(rls, 0.98, 0.25, 0.97)
lines(res3, col = "red")


######################################################
# An exmaple to show boundary correction.
# Load in two blaze-removed neighboring orders.
load("examples/BoundaryCorrection.Rdata")
plot(x1, type = "l", xlim=c(4610, 4710))
lines(x2, col="green")

corrected <- BC(x1, x2)
lines(corrected[[1]], col="cyan")
lines(corrected[[2]], col="blue")


######################################################
# An exmaple of continuous opacity correction is shown in 
# function/ContinuousOpacityCorrection.R


