# testing prior distributions



# want a shifted gamma distribution

# first use it as gamma to calculate. 

limits <- c(1.6, 3)
limits_tester <- limits-1

r0_gamma_pars <- get.gamma.par(p = c(0.025, 0.975), q = limits_tester, 
                               show.output = F, plot = F)

# x <- c(1:100)
# y <- rgamma(n=100, shape = r0_gamma_pars[1], rate = r0_gamma_pars[2])
# #plot(x, y)


# sucsceptibility, from Bens paper. 
# taking 95% to be lowest and median value in 2+ year olds (0.19, 0.5) <- immunity
# taking mean of other medians as our median = 0.3925
limits <- c(0.5,0.3925,0.19)
limits_sus <- 1 - limits
# fit to a beta
sus_beta_pars <- get.beta.par(p = c(0.025, 0.5, 0.975), 
                              q = limits_sus, 
                              show.output = F, 
                              plot = F)

#  x <- c(1:1000)
#  y <- rbeta(n=1000, shape1 = sus_beta_pars[1], shape2 = sus_beta_pars[2])
#  plot(x, y)
# hist(y, breaks=20)
