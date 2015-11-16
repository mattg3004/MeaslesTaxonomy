####################################

list <- structure(NA,class="result")
"[<-.result" <- function(x,...,value) {
    args <- as.list(match.call())
    args <- args[-c(1:2,length(args))]
    length(value) <- length(args)
    for(i in seq(along=args)) {
        a <- args[[i]]
        if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
    }
    x
}

########################################################################

##' Function to output animation ready data.
##' This allows the entry at each year use case data from future years

########################################################################
animation.gaussian.weighting.method <- function(window.length, regions, 
                                                  gaussian.st.dev, cutoff = 50, 
                                                  interp.resolution = 20){
    
    require(stats)
    list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = get.data.for.animation(regions)
    
    x = seq(1980, 2014) 
    
    ##' interpolate the datasets to have entries for all points in time once the interpolation is done.
    list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                        subset.vaccination, 
                                                                                                        subset.birth.rates, 
                                                                                                        subset.pop.by.year,
                                                                                                        x,
                                                                                                        x)
    
    ##' output matrices the correct size for our animation
    list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = prepare.matrices.for.animation(interp.subset.data, subset.data)
    
    ##' number of unique years that we will have data for. The longer the window, the less unique years of data.
    num.windows = 26 + (10 - window.length)
    
    ##' first year of data
    year = 1980
    
    ##' setting up the datasets
    coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
    incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
    mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
    mean.vac = matrix(0, length(subset.data[ , 1]), num.windows)
    
    ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
    ##' mean vaccination rate over periods of length given by the window length.
    for ( j in 1 : num.windows){
        for ( i in 1 : length(subset.data[ , 1])){
            coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
                mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
                if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
                    coeff.var[i, j]  =  0
                } 
            }
            incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                                 as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
            if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
                mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            }
            if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
                mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            }
        }
        year = year + 1
    }
    
    ##' we take the weighted average of the values that we have calculated for each year, where the weights are gaussian,
    ##' with a specified number of years for the standard deviation
    x1 = seq(1, length(coeff.var[1, ]))
    w1 = matrix(0, length(x1), length(x1))
    for (i in 1 : length(x1)){
        ##' set up the gaussian weights for averaging
        
        w.input = x1 - x1[i]
        w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
        
        ##' make sure that the weights add up to 1 for each of the specific weightings
        
        w1[i, ] = w1[i, ] / sum(w1[i, ])
    }
    
    ##' make a set of matrices that are the same size as the matrices containing the data.
    coeff.2 = coeff.var
    incidence.2 = incidence.per.1000
    mbr2 = mean.br
    mvacc2 = mean.vac
    for(i in 1 : length(coeff.var[1, ])){
        for(j in 1 : length(coeff.var[, 1])){
            ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
            coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
            incidence.2[j, i] = sum(incidence.per.1000[j, ] * w1[i, ], na.rm = T)
            mbr2[j, i] = mean.br[j, i] 
            mvacc2[j, i] = mean.vac[j, i]
            
            ##' Should we do weighted average of birth rate and vaccination rate?
            ##' If so uncomment the next two lines
            
            # mbr2[j, i] = sum(mean.br[j, ] * w1[i, ], na.rm = T)
            # mvacc2[j, i] = sum(mean.vac[j, i] * w1[i, ], na.rm = T)
        }
    }
    
    ##' set the original data to be equal to the weighted data
    coeff.var.cases = coeff.2
    incidence.per.1000 = incidence.2
    mean.br = mbr2
    mean.vac = mvacc2
    
    
    ##' set up the timeline on which we do the interpolation. 
    ##' The number of sections that the yearly data is split up to is given by interp.resolution  
    x = seq(1980 + (window.length - 1), 2014)
    xout = seq(1980 + (window.length - 1), 2014, 1/interp.resolution)
    
    ##' interpolate the data and add columns that contain the correspondin country and WHO region of each line
    
    coeff.var.cases = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(coeff.var.cases, x, xout))
    incidence.per.1000 = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(incidence.per.1000, x, xout))
    mean.br = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.br, x, xout))
    mean.vac = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.vac, x, xout))
    
    ##' round the data to 2 decimal places for ease of reading.
    mean.vac[, -(1:2)] = round(as.numeric(mean.vac[, -(1:2)]), 2)
    mean.br[, -(1:2)] = round(as.numeric(mean.br[, -(1:2)]), 2)
    incidence.per.1000[, -(1:2)] = round(as.numeric(incidence.per.1000[, -(1:2)]), 2)
    coeff.var.cases[, -(1:2)] = round(as.numeric(coeff.var.cases[, -(1:2)]), 2)
    
    
    ##' set up the output to be the appropriate size and add column labels.
    ##' Additionally add enough room to include additional data for each year that will be used 
    ##' to calibrate the data for each year, so that the minimum and maximum of vaccination rate is 0 and 100 each time.
    ##' This ensures that the colour scale is constant
    
    output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
    output.data  =  data.frame(output.data)
    colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
    
    ##' input the appropriate data to the outputs 
    output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 1], length(coeff.var.cases[1, -(1:2)]))
    output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 2], length(coeff.var.cases[1, -(1:2)]))
    count = 1
    for(i in 3 : length(coeff.var.cases[1, ])){
        output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
        output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
        output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
        output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
        output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
        count = count + 1
    }
    
    ##' add in the dummy data for each year to keep the scales constant.
    year.mins = matrix(0, length(xout), 2)
    
    for(i in 1 : length(xout)){
        t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
        year.mins[i, 1] = xout[i]
        year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
    }
    
    l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
    
    for(i in 1 : (2 * length(regions) * length(xout))){
        y = l[i, 6]
        j = which(year.mins[, 1] == y)
        l[i, 5]  =  year.mins[j, 2]
    }
    
    output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
    
    for( i in 1 : (2 * length(regions) * length(xout))){
        output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
        output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
    }
    
    ##' make sure that each column that should be numeric is numeric.
    output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
    output.data$Incidence  =  as.numeric(output.data$Incidence)
    output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
    output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
    output.data$Year   =  as.numeric(output.data$Year)
    output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
    return(output.data)
}




######################################################

##' Function to output animation ready data.
##' Output data for animation when we calculate incidence and cv for previous cases only

######################################################
animation.gaussian.weighting.method.back.only <- function(window.length, regions, 
                                                            gaussian.st.dev, cutoff = 50, 
                                                            interp.resolution = 20){
    
    require(stats)
    list[subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data] = get.data.for.animation(regions)
    
    x = seq(1980, 2014) 
    
    ##' interpolate the datasets to have entries for all points in time once the interpolation is done.
    list[interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop] = interp.datasets(subset.data, 
                                                                                                        subset.vaccination, 
                                                                                                        subset.birth.rates, 
                                                                                                        subset.pop.by.year,
                                                                                                        x,
                                                                                                        x)
    
    ##' output matrices the correct size for our animation
    list[mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac] = prepare.matrices.for.animation(interp.subset.data, subset.data)
    
    ##' number of unique years that we will have data for. The longer the window, the less unique years of data.
    num.windows = 26 + (10 - window.length)
    
    ##' first year of data
    year = 1980
    
    ##' setting up the datasets
    coeff.var = matrix(0, length(subset.data[ , 1]), num.windows)
    incidence.per.1000 = matrix(0, length(subset.data[ , 1]), num.windows)
    mean.br = matrix(0, length(subset.data[ , 1]), num.windows)
    mean.vac = matrix(0, length(subset.data[ , 1]), num.windows)
    
    ##' do calculations that calculate the coefficient of variation, incidence per 100, mean birth rate and
    ##' mean vaccination rate over periods of length given by the window length.
    for ( j in 1 : num.windows){
        for ( i in 1 : length(subset.data[ , 1])){
            coeff.var[i, j]  =  sd(interp.subset.data[i, paste(seq(year, year + window.length - 1))], na.rm = TRUE) /  
                mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            if(is.na( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)) == FALSE ){
                if( mean(as.numeric(interp.subset.data[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) == 0) {
                    coeff.var[i, j]  =  0
                } 
            }
            incidence.per.1000[i, j]  =  sum(as.numeric(interp.subset.data[i,  paste(seq(year, year + window.length - 1))]) / 
                                                 as.numeric(interp.subset.pop[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE) * 1000
            if(length(which(is.na(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))])))) < window.length){
                mean.br[i, j]  =  mean(as.numeric(interp.subset.br[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            }
            if(length(which(is.na(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))])))) < window.length){
                mean.vac[i, j]  =  mean(as.numeric(interp.subset.vacc[i, paste(seq(year, year + window.length - 1))]), na.rm = TRUE)
            }
        }
        year = year + 1
    }
    
    ##' we take the weighted average of the values that we have calculated for each year, where the weights are gaussian,
    ##' with a specified number of years for the standard deviation
    x1 = seq(1, length(coeff.var[1, ]))
    w1 = matrix(0, length(x1), length(x1))
    for (i in 1 : length(x1)){
        ##' set up the gaussian weights for averaging
        
        w.input = x1 - x1[i]
        w1[i, ] = output.weights.gaussian.with.cutoff(w.input, gaussian.st.dev, cutoff)
        
        ##' make sure that the weights add up to 1 for each of the specific weightings
        j = which(w.input == 0)
        w1[i,j : length(w1[i,])] = 0
        w1[i, ] = w1[i, ] / sum(w1[i, ])
    }
    
    ##' make a set of matrices that are the same size as the matrices containing the data.
    coeff.2 = coeff.var
    incidence.2 = incidence.per.1000
    mbr2 = mean.br
    mvacc2 = mean.vac
    for(i in 1 : length(coeff.var[1, ])){
        for(j in 1 : length(coeff.var[, 1])){
            ##' make the entries of these newly created matrices to be the weighted averages of the originally calculated datasets
            coeff.2[j, i] = sum(coeff.var[j, ] * w1[i, ], na.rm = T)
            incidence.2[j, i] = sum(incidence.per.1000[j, ] * w1[i, ], na.rm = T)
            mbr2[j, i] = mean.br[j, i] 
            mvacc2[j, i] = mean.vac[j, i]
            
            ##' Should we do weighted average of birth rate and vaccination rate?
            ##' If so uncomment the next two lines
            
            # mbr2[j, i] = sum(mean.br[j, ] * w1[i, ], na.rm = T)
            # mvacc2[j, i] = sum(mean.vac[j, i] * w1[i, ], na.rm = T)
        }
    }
    
    ##' set the original data to be equal to the weighted data
    coeff.var.cases = coeff.2
    incidence.per.1000 = incidence.2
    mean.br = mbr2
    mean.vac = mvacc2
    
    
    ##' set up the timeline on which we do the interpolation. 
    ##' The number of sections that the yearly data is split up to is given by interp.resolution  
    x = seq(1980 + (window.length - 1), 2014)
    xout = seq(1980 + (window.length - 1), 2014, 1/interp.resolution)
    
    ##' interpolate the data and add columns that contain the correspondin country and WHO region of each line
    
    coeff.var.cases = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(coeff.var.cases, x, xout))
    incidence.per.1000 = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(incidence.per.1000, x, xout))
    mean.br = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.br, x, xout))
    mean.vac = cbind(interp.subset.data[,(1:2)], interpolate.give.dataset(mean.vac, x, xout))
    
    ##' round the data to 2 decimal places for ease of reading.
    mean.vac[, -(1:2)] = round(as.numeric(mean.vac[, -(1:2)]), 2)
    mean.br[, -(1:2)] = round(as.numeric(mean.br[, -(1:2)]), 2)
    incidence.per.1000[, -(1:2)] = round(as.numeric(incidence.per.1000[, -(1:2)]), 2)
    coeff.var.cases[, -(1:2)] = round(as.numeric(coeff.var.cases[, -(1:2)]), 2)
    
    
    ##' set up the output to be the appropriate size and add column labels.
    ##' Additionally add enough room to include additional data for each year that will be used 
    ##' to calibrate the data for each year, so that the minimum and maximum of vaccination rate is 0 and 100 each time.
    ##' This ensures that the colour scale is constant
    
    output.data = matrix(0, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + (2 * length(regions) * length(xout)), 7)
    output.data  =  data.frame(output.data)
    colnames(output.data) = c("Country", "Coefficient.of.Variation", "Incidence", "Mean.vaccination", "Mean.birth.rate", "Year", "WHO_REGION")
    
    ##' input the appropriate data to the outputs 
    output.data$Country[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 1], length(coeff.var.cases[1, -(1:2)]))
    output.data$WHO_REGION[seq(1, (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]))] = rep(interp.subset.data[, 2], length(coeff.var.cases[1, -(1:2)]))
    count = 1
    for(i in 3 : length(coeff.var.cases[1, ])){
        output.data$Coefficient.of.Variation[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = coeff.var.cases[, i]
        output.data$Incidence[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = incidence.per.1000[, i]
        output.data$Mean.vaccination[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.vac[, i]
        output.data$Mean.birth.rate[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = mean.br[, i]
        output.data$Year[seq(((count - 1) * (length(coeff.var.cases[, 1]))) + 1, count * (length(coeff.var.cases[, 1])))] = xout[i - 2]  
        count = count + 1
    }
    
    ##' add in the dummy data for each year to keep the scales constant.
    year.mins = matrix(0, length(xout), 2)
    
    for(i in 1 : length(xout)){
        t  =  subset(output.data, output.data$Year ==  unique(xout)[i])
        year.mins[i, 1] = xout[i]
        year.mins[i, 2] = as.numeric(min(t$Mean.birth.rate,na.rm = T)  )
    }
    
    l =  expand.grid("", -1, 0, c(0,100), 0, xout, regions)
    
    for(i in 1 : (2 * length(regions) * length(xout))){
        y = l[i, 6]
        j = which(year.mins[, 1] == y)
        l[i, 5]  =  year.mins[j, 2]
    }
    
    output.data[((length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + 1 ): length(output.data[, 1]), ]  =  l
    
    for( i in 1 : (2 * length(regions) * length(xout))){
        output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] = regions[as.numeric(output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)]) + i, 7] )]
        output.data[ (length(coeff.var.cases[, 1])) * length(coeff.var.cases[1, -(1:2)])  + i, 1] = ""
    }
    
    ##' make sure that each column that should be numeric is numeric.
    output.data$Coefficient.of.Variation = as.numeric(output.data$Coefficient.of.Variation)
    output.data$Incidence  =  as.numeric(output.data$Incidence)
    output.data$Mean.vaccination  =  as.numeric(output.data$Mean.vaccination)
    output.data$Mean.birth.rate  =  as.numeric(output.data$Mean.birth.rate)
    output.data$Year   =  as.numeric(output.data$Year)
    output.data$Coefficient.of.Variation[which(output.data$Coefficient.of.Variation == "Inf")] = 0
    return(output.data)
}



############################

##' Import data needed to generate animations

############################
get.data.for.animation <- function(regions){
    
    cases.by.country.by.year = read.csv("data/Measles_cases_by_year.csv", stringsAsFactors = FALSE)
    #cases.by.country.by.year = read.csv("Measles_cases_by_year2.csv", stringsAsFactors = FALSE)
    Birth.rates = read.csv("data/Birth_rates.csv", stringsAsFactors = FALSE)
    pop.by.year = read.csv("data/All_populations.csv", stringsAsFactors = FALSE)
    vacc.rates = read.csv("data/Measles_vac_all.csv", stringsAsFactors = FALSE)
    
    Birth.rates$X2013 = Birth.rates$X2012
    Birth.rates$X2014 = Birth.rates$X2013
    pop.by.year$X2014 = pop.by.year$X2013
    vacc.rates$X2014 = vacc.rates$X2013
    
    subset.pop.by.year = subset(pop.by.year, pop.by.year$WHO_REGION %in% regions)
    subset.vaccination = subset(vacc.rates, vacc.rates$WHO_REGION %in% regions)
    subset.birth.rates = subset(Birth.rates, Birth.rates$WHO_REGION %in% regions)
    subset.data = subset(cases.by.country.by.year, cases.by.country.by.year$WHO_REGION %in% regions)
    
    missing1 = setdiff(subset.vaccination$Country, subset.pop.by.year$Country.Name)
    if(length(missing1) > 0){
        j = which(subset.vaccination$Country %in% missing1)
        subset.vaccination = subset.vaccination[-j, ]
    }
    missing2 = setdiff(subset.birth.rates$Country, subset.vaccination$Country)
    if(length(missing2) > 0){
        j = which(subset.birth.rates$Country %in% missing2)
        subset.birth.rates = subset.birth.rates[-j, ]
    }
    missing3 = setdiff(subset.data$Cname, subset.vaccination$Country)
    if(length(missing3) > 0){
        j = which(subset.data$Cname %in% missing3)
        subset.data = subset.data[-j, ]
    }
    
    
    p1  =  subset.pop.by.year
    p2  =  subset.vaccination
    p3  =  subset.birth.rates
    p4  =  subset.data
    
    for ( i in 1 : length(subset.vaccination[, 1])){
        C  =  subset.vaccination$Country[i]
        p2[i, ]  =  subset.vaccination[i, ]
        j = which(subset.pop.by.year$Country.Name == C)
        p1[i, ]  =  subset.pop.by.year[j, ]
        j = which(subset.birth.rates$Country == C)
        p3[i, ]  =  subset.birth.rates[j, ]
        j = which(subset.data$Cname == C)
        p4[i, ]  =  subset.data[j, ]
    }
    subset.pop.by.year = p1
    subset.vaccination = p2
    subset.birth.rates = p3
    subset.data = p4
    
    return(list(subset.pop.by.year, subset.vaccination, subset.birth.rates, subset.data))
}




############################

##' Get datasets into the same format

############################

interp.datasets <- function(subset.data, 
                            subset.vaccination, 
                            subset.birth.rates, 
                            subset.pop.by.year, 
                            x,
                            xout){
    
    interp.subset.data = matrix(0, length(subset.data[, 1]), length(xout) + 2)
    interp.subset.data[, 1] = subset.data$Cname
    interp.subset.data[, 2] = subset.data$WHO_REGION
    
    interp.subset.vacc = matrix(0, length(subset.vaccination[, 1]), length(xout) + 2)
    interp.subset.vacc[, 1] = subset.data$Cname
    interp.subset.vacc[, 2] = subset.data$WHO_REGION
    
    interp.subset.br = matrix(0, length(subset.birth.rates[, 1]), length(xout) + 2)
    interp.subset.br[, 1] = subset.data$Cname
    interp.subset.br[, 2] = subset.data$WHO_REGION
    
    interp.subset.pop = matrix(0, length(subset.pop.by.year[, 1]), length(xout) + 2)
    interp.subset.pop[, 1] = subset.data$Cname
    interp.subset.pop[, 2] = subset.data$WHO_REGION
    
    for ( i in 1 : length(subset.pop.by.year[, 1])){
        y = subset.data[i, paste("X", seq(1980, 2014), sep = "")]
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
        interp.subset.data[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
        colnames(interp.subset.data) = c("Country", "WHO_REGION", seq(1980, 2014))
        
        y1 = subset.vaccination[i, paste("X", seq(1980, 2014), sep = "")]
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y1),  method = "linear", xout )
        interp.subset.vacc[i, 3: length(interp.subset.vacc[1, ])]  =  round(ww, 2)
        colnames(interp.subset.vacc) = c("Country", "WHO_REGION", seq(1980, 2014))
        
        y2 = subset.birth.rates[i, paste("X", seq(1980, 2014), sep = "")]
        if(length(which(is.na(y2) == F)) < 2) 
        {interp.subset.br[i, 3: length(interp.subset.data[1, ])]  = 0} else{
            list[qq,ww] =  approx (as.numeric(x), as.numeric(y2),  method = "linear", xout )
            interp.subset.br[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
        }
        colnames(interp.subset.br) = c("Country", "WHO_REGION", seq(1980, 2014))
        
        
        y3 = subset.pop.by.year[i, paste("X", seq(1980, 2014), sep = "")]
        list[qq,ww] =  approx (as.numeric(x), as.numeric(y3),  method = "linear", xout )
        interp.subset.pop[i, 3: length(interp.subset.data[1, ])]  =  round(ww, 2)
        colnames(interp.subset.pop) = c("Country", "WHO_REGION", seq(1980, 2014))
    }
    
    return(list(interp.subset.data, interp.subset.vacc, interp.subset.br, interp.subset.pop))
    
}




############################

##' Continue preparing matrices

############################

prepare.matrices.for.animation <- function(interp.subset.data, subset.data){
    
    mean.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
    mean.cases[, 1] = subset.data$Cname
    mean.cases[, 2] = subset.data$WHO_REGION
    
    coeff.var.cases = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
    coeff.var.cases[, 1] = subset.data$Cname
    coeff.var.cases[, 2] = subset.data$WHO_REGION
    
    incidence.per.1000 = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
    incidence.per.1000[, 1] = subset.data$Cname
    incidence.per.1000[, 2] = subset.data$WHO_REGION
    
    mean.br = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
    mean.br[, 1] = subset.data$Cname
    mean.br[, 2] = subset.data$WHO_REGION
    
    mean.vac = matrix(0, length(interp.subset.data[, 1]), length(interp.subset.data[1, ]))
    mean.vac[, 1] = subset.data$Cname
    mean.vac[, 2] = subset.data$WHO_REGION
    
    return(list(mean.cases, coeff.var.cases, incidence.per.1000, mean.br, mean.vac))
}





############################

##' output the gaussian weight for the averaging process

############################
output.weights.gaussian.with.cutoff <- function(x, st.dev, cutoff, neg.only = F){
    
    weights = dnorm(x, mean = 0, sd = st.dev)
    if(neg.only == T){
        weights[which(x > cutoff)] = 0
    } else{
        weights[which(abs(x) > cutoff)] = 0
    }
    return(weights)
}



############################

##' Interpolate datasets to the specified level

############################

interpolate.give.dataset <- function(data,  
                                     x,
                                     xout){
    interp.data = matrix(0, length(data[, 1]), length(xout))
    for ( i in 1 : length(data[, 1])){
        y = data[i, ]
        if(length(which(!is.na(y)) == F) < 2){
            interp.data[i, ]  =  0} else{
                list[qq,ww] =  approx (as.numeric(x), as.numeric(y),  method = "linear", xout )
                #interp.data[i, ]  =  round(ww, 2)
                interp.data[i, ]  =  ww
            }
    }
    
    return(interp.data)
    
}



############################

##' Output a dendrogram by specifiying the data to use along with the WHO regions and 
##' the variables from the dataset to use for the distance calculation

############################


make.dendrogram <- function(data, regions, variables){
    require(ape)
    
    ##' Get a list of the countries in the specified regions, and discount any with empty names
    Countries = as.character(unique(data$Country[which(data$WHO_REGION %in% regions)]))
    Countries = Countries[which(Countries != "")]
    
    
    ##' discount first year in the dataset and look at only yearly values, rather than all interpolated time points
    years = seq(min(data$Year) + 1, max(data$Year))
    Z = data[which(data$Year %in% years), ]
    
    
    ##' set up matrices to contain data specified by variables input
    
    data.matrix = matrix(0, length(Countries), length(years) * length(variables))
    data.list = matrix(0, length(Countries) * length(years), length(variables))
    rownames(data.matrix) = Countries
    
    
    ##' Put data in the data.matrix
    ##' Also keep a list for each variable, so that we can later calculate
    ##' the mean and variance for the separate variables
    count1 = 1
    count2 = 1
    
    for(k in 1 : length(variables)){
        for(j in 1 : length(years)){
            y = years[j]
            
            for(i in 1 : length(Countries)){
                c = Countries[i]
                data.matrix[i, count1] = Z[which(Z$Country == c & 
                                                     Z$Year == y), variables[k]]
                data.list[count2, k] = data.matrix[i, count1]
                count2 = count2 + 1
            }
            
            count1 = count1 + 1
        }
        count2 = 1
    }
    
    
    ##' normalise the data
    
    for(i in 1 : length(variables)){
        A = data.matrix[, (((i-1)*length(years)) + 1):(i*length(years))]
        B = data.list[, i]
        data.matrix[, (((i-1)*length(years)) + 1):(i*length(years))] = (A - mean(B))/as.numeric(sqrt(var(B)))
    }
    ##' join incidence and cv data together. Then calculate the distance between them and plot a dendrogram
    
    hc = hclust(dist(data.matrix))
    #plot(hc)
    
    op = par(bg = "#DDE3CA")
    #op = par(bg = rgb(1,1,1))
    
    plot(as.phylo(hc), cex = 0.7, label.offset = 0.1)
}
