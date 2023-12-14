The workflow for the R Markdown scripts is from 00 to 04:

00-isocape-data-preparation.Rmd
01-isocape-model-fits-regressions.Rmd
02-isocape-kriging-sinusoidal-parameters.Rmd
03-isocape-monthly-residuals.Rmd
04-isocape-jack-knife.Rmd

Associated with each of these is an html file with the R code, output, and documentation.
Except for 04 which takes a very long time to run (a week or more).

Output data at each step is saved, and feeds into the following steps. 


Data Files
----------

Raw input data (e.g. isotope values at sites, site information) are in the directory "data"

Some input data that required substantial pre-processing time (mostly the VCSN climate data) 
are in the directory "ProcessedData"

Output data at each step is saved, and feeds into the following steps. These intermediate
data are saved as R objects in "Output/Data"

Results at each step are saved in sub-directories under "Output":

ModelFitResults   (results from 01)
KrigingResults    (results from 02)
MonthlyResults    (results from 03)
JackKnifeResults  (results from 04)





