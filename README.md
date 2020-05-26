# Multiple imputation of missing data
R codes for statistical and econometrics models applied to transport data

This code has been produced as a part of my post-doctoral research. Please cite the following article if you use this code in any kind:

Afghari, A.P., Washington, S., Prato, C. and Haque, M.M., 2019. Contrasting case-wise deletion with multiple imputation and latent variable approaches to dealing with missing observations in count regression models. Analytic Methods in Accident Research, 24, p.100104.

This R code corresponds with multiple imputation of missing data in regression models. The built-in function for multiple imputation in R applies the multiple imputation on data but stops after imputing missing values multiple times. This is simply repeating the imputation multiple times. The definition of 'multiple imputation' in the statistical literature goes beyond naive repetition of imputation multiple times; it rather provides an estimate of the modelling uncertainty created by the data ‘missing-ness’, as distinct from the natural variation in the data. This code provides the latter information on the variability of estimates caused by missing data in regression models.

Although this code has been written for multiple imputation in count regression models, it may be adjusted for any types of regression by replacing the likelihood function with the likelihood function of the preferred model.

Prior to running the code, make sure to have "maxLik", "randtoolbox", and "MICE" packages installed.
