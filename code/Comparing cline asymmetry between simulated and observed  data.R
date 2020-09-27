# read data file containing position of mating trait and neutral cline centers by the end (generation 500) of 8 simulated scenario
# Only clines that moved towards dominant homozygote are presented
sim.cline.param <- read.table("simulated.cline.parameters.txt", header=T)

# assign dispersal distance value as was used in simulations
meandispersal <- 0.01

# assign empty vector
vec <- c()

# loop to generate 10,000,000 random dispersal events
for(i in 1:10000000){
  vec[i] <- rnorm(1, mean=0, sd=meandispersal)
}

# calculate absolute mean value of simulated dispersal distance
empirical_mean_abs <- mean(abs(vec))

# calculate desplacement between mating and neutral clines across 8 simulated scenario
sim.cline.param$center.displacement <- sim.cline.param$Mating.center - sim.cline.param$Neutral.center

# convertion factor to wagtail units is based on average post natal dispersal distance of Motacilla alba yarrellii of 20.34km (Dougall 1991).
# Dougall, T. W. (1991). Winter distribution and associated movements of northern Pied Wagtails Motacilla alba yarrellii, as shown by ringing. Ringing & Migration, 12, 1–15.
conversion_factor <- empirical_mean_abs / 20.34 

# convert values of simulated cline displacement to wagtail units
sim.cline.param$center.displacement.in.wagtail.units <- sim.cline.param$center.displacement/conversion_factor

mean(sim.cline.param$center.displacement.in.wagtail.units)
sd(sim.cline.param$center.displacement.in.wagtail.units)
