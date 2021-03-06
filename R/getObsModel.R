#' Create the observation model component of a sparta JAGS model 
#' 
#' This function is primarily for internal use within \code{getModelFile}. It is used to 
#' write an observation model that fits the users needs. The model is returned as a character.
#' 
#' @param modeltype Character, one of: jul_date, catlistlength, contlistlength.
#' See \code{occDetFunc} for more information.
#' @param verbose Logical, if true progress is reported to the console
#' @return A character, of JAGS model code, that describes the observation model.
#' @export

getObsModel <- function(modeltype, verbose = FALSE){
  
  basemodel <- "logit(p[j]) <-  alpha.p[Year[j]]"
  addVar <- NULL
  
  # if list length (categorical) is not a specification then we 
  
  for(mtype in modeltype){
    
    switch(EXPR = tolower(mtype),
           
           jul_date = {
             basemodel <- paste(
               basemodel,
               'beta3*f_JD[JulDate[j]]',
               sep = ' + ')
             addVar <- c(addVar,
                         "beta1 ~ dunif(1, 366)\nbeta2 ~ dt(0, 1, 1)T(0,)\nbeta3 ~ dnorm(0, 0.0001)\nfor (d in 1:366){\nf_JD[d] <- 1/((2*3.14159265359)^0.5 * beta2) * exp(-((d - beta1)^2 / (2* beta2^2)))\n}\n")
           },
           
           catlistlength = {
             basemodel <- paste(
               basemodel,
               'dtype2.p*DATATYPE2[j] + dtype3.p*DATATYPE3[j]',
               sep = ' + ')
             addVar <- c(addVar,
                         "dtype2.p ~ dnorm(0, 0.01)\ndtype3.p ~ dnorm(0, 0.01)\n") 
           },
           
           contlistlength = {
             basemodel <- paste(
               basemodel,
               'LL.p*logL[j]',
               sep = ' + ')
             addVar <- c(addVar,
                         "LL.p ~ dunif(dtype2p_min, dtype2p_max)\n")
             },
           
           mixlistlength = {
             basemodel <- paste(
               basemodel,
               'LL.p*logL[j]',
               'dtype2.p*DATATYPE2[j] + dtype3.p*DATATYPE3[j]',
               sep = ' + ')
             addVar <- c(addVar,
                         "LL.p ~ dunif(dtype2p_min, dtype2p_max)\n",
                         "dtype2.p ~ dnorm(0, 0.01)\ndtype3.p ~ dnorm(0, 0.01)\n")
             },
           
           mixlistlength2 = {
             basemodel <- paste(
               "logit(p[j]) <-  alpha.p[Year[j]]*(1-DATATYPE3[j])",
               'LL.p*logL[j]',
               'dtype2.p*DATATYPE2[j] + dtype3.p*DATATYPE3[j]',
               sep = ' + ')
             addVar <- c(addVar,
                         "LL.p ~ dunif(dtype2p_min, dtype2p_max)\n",
                         "dtype2.p ~ dnorm(0, 0.01)\ndtype3.p ~ dnorm(0, 0.01)\n")
           }
    )
  }
  
  fullModel <- paste0(paste(addVar, collapse = '\n'),
                     paste0('### Observation Model\n',
                            'for(j in 1:nvisit) {\n',
                            '  y[j] ~ dbern(Py[j])\n', 
                            '  Py[j]<- z[Site[j],Year[j]]*p[j]'),
                     '\n  ', basemodel, '\n',
                     '}')
  
  if(verbose) cat(paste('Observation model:\n', fullModel, '\n'))
  
  return(fullModel)
  
}