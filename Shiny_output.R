source("Functions.R")
id.names = LETTERS[1:6]
regions = c("EMR","EUR","AFR","AMR","WPR","SEAR")

interp.resolution = 20
official.vacc = 0


gaussian.weighted.data = animation.gaussian.weighting.method.2(window.length = 4, regions,
                                                               gaussian.st.dev = 4, cutoff = 50, 
                                                               interp.resolution = 20)

anim.data = gaussian.weighted.data

gaussian.weighted.data.back.only = animation.gaussian.weighting.method.back.only(window.length = 10, regions,
                                                                                   gaussian.st.dev = 3, cutoff = 50, 
                                                                                   interp.resolution = 20)
anim.data = gaussian.weighted.data.back.only

shiny::runApp('shiny-plots-all-regions-weighted-means')


