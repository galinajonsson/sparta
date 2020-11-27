context("Test formatOccData")

# Create data
n <- 15000 #size of dataset
nyr <- 20 # number of years in data
nSamples <- 100 # set number of dates
nSites <- 50 # set number of sites
set.seed(125)

# Create somes dates
first <- as.Date(strptime("2010/01/01", "%Y/%m/%d")) 
last <- as.Date(strptime(paste(2010+(nyr-1),"/12/31", sep=''), "%Y/%m/%d")) 
dt <- last-first 
rDates <- first + (runif(nSamples)*dt)

# taxa are set as random letters
taxa <- sample(letters, size = n, TRUE)

# three sites are visited randomly
site <- sample(paste('A', 1:nSites, sep=''), size = n, TRUE)

# the date of visit is selected at random from those created earlier
survey <- sample(rDates, size = n, TRUE)

# set the closure period to be in 2 year bins
closure_period <- ceiling((as.numeric(format(survey,'%Y')) - 2009)/2)

# create survey variable that is not a date and a closure period to match
survey_numbered <- as.integer(as.factor(survey))


replicate <- rep(1, n)

test_that("Test formatOccData", {

  expect_warning(visitData <- formatOccData(taxa = taxa, site = site, survey = survey),
                 '871 out of 15000 observations will be removed as duplicates')
  
  head_spp_vis <- structure(list(visit = c("A102010-04-141", "A102010-04-221", "A102010-08-291", 
                            "A102010-11-041", "A102011-02-091", "A102011-03-091"), a = c(FALSE, 
                            FALSE, FALSE, FALSE, FALSE, FALSE), b = c(FALSE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE), c = c(FALSE, FALSE, FALSE, FALSE, FALSE, 
                            FALSE), d = c(FALSE, TRUE, FALSE, TRUE, FALSE, FALSE), e = c(FALSE, 
                            TRUE, FALSE, FALSE, FALSE, FALSE), f = c(FALSE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE), g = c(FALSE, FALSE, TRUE, FALSE, FALSE, 
                            FALSE), h = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), i = c(FALSE, 
                            FALSE, FALSE, FALSE, FALSE, FALSE), j = c(FALSE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE), k = c(FALSE, FALSE, FALSE, FALSE, FALSE, 
                            FALSE), l = c(TRUE, FALSE, FALSE, FALSE, FALSE, TRUE), m = c(FALSE, 
                            FALSE, FALSE, FALSE, FALSE, FALSE), n = c(FALSE, TRUE, TRUE, 
                            FALSE, FALSE, FALSE), o = c(TRUE, FALSE, FALSE, FALSE, FALSE, 
                            FALSE), p = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE), q = c(FALSE, 
                            TRUE, FALSE, FALSE, FALSE, FALSE), r = c(TRUE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE), s = c(FALSE, FALSE, FALSE, TRUE, TRUE, 
                            FALSE), t = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE), u = c(FALSE, 
                            FALSE, FALSE, FALSE, FALSE, FALSE), v = c(FALSE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE), w = c(TRUE, FALSE, FALSE, FALSE, FALSE, 
                            FALSE), x = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE), y = c(FALSE, 
                            FALSE, FALSE, FALSE, FALSE, FALSE), z = c(TRUE, FALSE, FALSE, 
                            FALSE, FALSE, FALSE)), .Names = c("visit", "a", "b", "c", "d", 
                            "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", 
                            "r", "s", "t", "u", "v", "w", "x", "y", "z"), row.names = c(NA, 
                            6L), class = "data.frame")
  

 head_occDetdata <- structure(list(visit = c("A102010-04-141", "A102010-04-221", 
                                             "A102010-08-291", "A102010-11-041",
                                             "A102011-02-091", "A102011-03-091"),
                                   site = c("A10", "A10", "A10", "A10", "A10",
                                            "A10"),
                                   L = c(5L, 5L, 2L, 2L, 1L, 2L),
                                   TP = c(2010, 2010, 2010, 2010, 2011, 2011)),
                              row.names = c(1L, 6L, 11L, 13L, 15L, 16L),
                              class = "data.frame")
 
  expect_identical(head(visitData$spp_vis), head_spp_vis)
  expect_identical(head(visitData$occDetdata), head_occDetdata)
    
})

test_that("Test formatOccData errors", {
  
  expect_error(visitData <- formatOccData(taxa = head(taxa), site = site, survey = survey, closure_period = closure_period, replicate = replicate),
               'The following arguements are not of equal length: taxa, site, survey, closure_period, replicate')
  expect_error(visitData <- formatOccData(taxa = taxa, site = head(site), survey = survey, closure_period = closure_period, replicate = replicate),
               'The following arguements are not of equal length: taxa, site, survey, closure_period, replicate')
  expect_error(visitData <- formatOccData(taxa = taxa, site = site, survey = head(survey), closure_period = closure_period, replicate = replicate),
               'The following arguements are not of equal length: taxa, site, survey, closure_period, replicate')
  expect_error(visitData <- formatOccData(taxa = taxa, site = site, survey = survey, closure_period=head(closure_period), replicate = replicate),
               'The following arguements are not of equal length: taxa, site, survey, closure_period, replicate')
  expect_error(visitData <- formatOccData(taxa = taxa, site = site, survey = survey, closure_period=closure_period, replicate = head(replicate)),
               'The following arguements are not of equal length: taxa, site, survey, closure_period, replicate')
  
})

test_that("Test formatOccData date requirement errors", {
  
  expect_warning(visitData <- formatOccData(taxa = taxa, site = site, survey = survey, closure_period = closure_period),
                 '871 out of 15000 observations will be removed as duplicates')
  expect_error(suppressWarnings(visitData <- formatOccData(taxa = taxa, site = site, survey = survey_numbered)),
               'survey must be a date if closure_period not supplied')
  expect_error(suppressWarnings(visitData <- formatOccData(taxa =taxa, site = site, survey = survey_numbered, includeJDay = TRUE, closure_period = closure_period)),
               'survey must be a date if Julian Date is to be included')
})

test_that("Test formatOccData specified closure period", {
  
  expect_warning(visitData <- formatOccData(taxa = taxa, site = site, survey = survey, closure_period = closure_period),
                 '871 out of 15000 observations will be removed as duplicates')
  
  head_occDetdata_cp <- structure(list(visit = c("A102010-04-141", "A102010-04-221", 
                                                 "A102010-08-291", "A102010-11-041",
                                                 "A102011-02-091", "A102011-03-091"),
                                       site = c("A10", "A10", "A10", "A10", "A10", "A10"),
                                       L = c(5L, 5L, 2L, 2L, 1L, 2L),
                                       TP = c(1L, 1L, 1L, 1L, 1L, 1L)),
                                  row.names = c(1L, 6L, 11L, 13L, 15L, 16L),
                                  class = "data.frame")
  
  expect_identical(head(visitData$occDetdata), head_occDetdata_cp)
  
})