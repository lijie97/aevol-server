library("ggplot2")
library("cowplot")

setwd("/Users/charlesrocabert/git/aevol/examples/wt/stats")

b = read.csv("stat_fitness_best.txt", skip=22, sep=" ", h=F)
m = read.csv("stat_fitness_glob.txt", skip=22, sep=" ", h=F)

b = b[1:(length(b[,1])-1),]
m = m[1:(length(m[,1])-1),]

names(b) = c("g", "popsize", "fitness", "log_sigma", "log_sigma_sd", "est_sigma", "est_sigma_sd", "genome_size", "error", "1", "2", "3", "4", "5", "6")
names(m) = c("g", "popsize", "fitness", "log_sigma", "log_sigma_sd", "est_sigma", "est_sigma_sd", "genome_size", "error", "1", "2", "3", "4", "5", "6")

SIGMA_MAX = 0.1

mean(b$fitness)

# p1 = ggplot() +
#   geom_line(data=m, aes(x=g, y=log10(fitness))) +
#   geom_line(data=b, aes(x=g, y=log10(fitness)), color="cornflowerblue") +
#   ggtitle("Fitness (log10)")
# p2 = ggplot() +
#   geom_line(data=m, aes(x=g, y=genome_size)) +
#   geom_line(data=b, aes(x=g, y=genome_size), color="cornflowerblue") +
#   ggtitle("Genome size")
# p3 = ggplot() +
#   geom_line(data=m, aes(x=g, y=log_sigma)) +
#   geom_line(data=b, aes(x=g, y=log_sigma), color="cornflowerblue") +
#   ggtitle("Log-sigma")
# p4 = ggplot() +
#   geom_line(data=m, aes(x=g, y=log_sigma_sd)) +
#   geom_line(data=b, aes(x=g, y=log_sigma_sd), color="cornflowerblue") +
#   ggtitle("Log-sigma variability")
# p5 = ggplot() +
#   geom_line(data=m, aes(x=g, y=est_sigma)) +
#   geom_line(data=b, aes(x=g, y=est_sigma), color="cornflowerblue") +
#   ggtitle("Estimated sigma")
# p6 = ggplot() +
#   geom_line(data=m, aes(x=g, y=est_sigma_sd)) +
#   geom_line(data=b, aes(x=g, y=est_sigma_sd), color="cornflowerblue") +
#   ggtitle("Estimated sigma variability")
# 
# plot_grid(p1, p2, p3, p4, p5, p6, labels="AUTO", ncol=2)

