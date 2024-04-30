# Functions for simulating metabolism data
# A Carter
# 4/2022

# litterfall ####
calc_litter_oneyear <- function(annual_litter, LAI){
    LAI <- zoo::na.approx(LAI, na.rm = F)
    diffLAI <- c(0, diff(LAI))
    diffLAI <- case_when(diffLAI < 0 ~ diffLAI,
                         TRUE ~ 0)
    litter <- annual_litter * diffLAI/sum(diffLAI, na.rm = T)

    return(litter)
}

calc_litter_from_LAI <- function(dd){
    if(!all(c('date', 'LAI') %in% colnames(dd))){
        print('Your dataframe must have columns called \'date\' and \'LAI\'')
        return()
    }
    if(!'year' %in% colnames(dd)){
        dd$year = lubridate::year(dd$date)
    }
    tmp <- data.frame()
    for(y in unique(dd$year)){
        yy <- dd[dd$year == y,]
        yy$litter = calc_litter_oneyear(100*max(yy$LAI), yy$LAI)
        tmp <- bind_rows(tmp, yy)
    }

    return(tmp$litter)
}

# Antecedent GPP ####
calc_antecedent_GPP <- function(GPP, nweights, w = NULL){
    GPP <- zoo::na.approx(GPP, na.rm = F)
    if(sum(is.na(GPP)) > 0) {
        print('GPP should not start or end with NA')
        return()
    }

    if(is.null(w)) w = rep(1/nweights, nweights)
    Pant <- GPP <- c(rep(GPP[1], nweights), GPP)
    ndays <- length(GPP)

    for(i in (nweights+1):ndays){
        Pvec <- numeric(nweights)
        for(j in 1:nweights){
            Pvec[j] <- w[j]*GPP[i-j]
        }
        Pant[i]<-sum(Pvec)
    }

    Pant <- Pant[(nweights+1):ndays]

    return(Pant)
}

calc_antecedent_GPP_vi <- function(GPP, nweights, w = NULL,
                                   intervals = NULL, scale = NULL){

    GPP <- zoo::na.approx(GPP, na.rm = F)
    if(sum(is.na(GPP)) > 0) {
        print('GPP should not start or end with NA')
        return()
    }

    if(is.null(w)) w = rep(1/nweights, nweights)
    if(is.null(intervals)){
        index <- seq(0,(nweights-1), by = 1)
        intervals = scale * 2 ^ index
    }

    n_pre <- sum(intervals)
    Pant <- GPP <- c(rep(GPP[1], n_pre), GPP)
    ndays <- length(GPP)

    for (i in (n_pre+1):ndays){
        Pvec <- numeric(nweights)
        N <- 0
        for(j in 1:nweights){
            n <- intervals[j]
            pp <- numeric(n)
            for(k in 1:n){
                pp[k] <- GPP[i-N-k]
            }
            Pvec[j] <- w[j] * mean(pp)
            N <- N + n
        }
        Pant[i]<-sum(Pvec)
    }

    Pant <- Pant[(n_pre+1):ndays]

    return(Pant)
}

# Simulate detrital carbon and detrital based respiration ####
calc_rate_coef <- function(temp_C,
                           K_20 = 0.01,# respiration rate coefficient at 20C
                           E_a = 0.63 # activation energy for heterotrophic respiration (eV)
                           ){
    k_b = 8.6173*10^-5 # boltzmann's constant in eV/K
    K = K_20 * exp(-E_a/k_b *(1/(temp_C+273) - 1/293))

    return(K)
}

