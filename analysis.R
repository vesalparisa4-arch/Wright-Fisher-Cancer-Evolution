# Extracted from original Tutorial 4 script, part 1.
if (!requireNamespace('ggplot2', quietly=TRUE)) install.packages('ggplot2', repos='https://cloud.r-project.org')
library(ggplot2)
out <- 'Week4_R_outputs'
dir.create(out,showWarnings=FALSE)

# 1. Wright-Fisher process (neutral, fixed n, no mutation, no selection)
wright_fisher <- function(n, ngen=100, nrep=50, gamma=0.1) {
  stopifnot(n>=1, ngen>=0, nrep>=1, gamma>=0, gamma<=1)
  # Start each replicate at the integer cell count closest to n * gamma.
  init <- round(n * gamma)
  states <- matrix(NA_real_, nrow=ngen+1, ncol=nrep)
  states[1,] <- init/n
  if (ngen>0) for (g in seq_len(ngen)) {
    states[g+1,] <- rbinom(nrep, size=n, prob=states[g,])/n
  }
  data.frame(n=n, gamma=gamma, generation=rep(0:ngen, times=nrep),
             replicate=rep(seq_len(nrep), each=ngen+1), proportion=as.vector(states))
}
set.seed(20261123)
params <- expand.grid(n=c(50,100,1000,5000), gamma=c(.01,.1,.5,.8))
sim_list <- lapply(seq_len(nrow(params)), function(i) {
  wright_fisher(params$n[i], 100, 50, params$gamma[i])
})
sim <- do.call(rbind,sim_list)
write.csv(sim, file.path(out,'01_WF_simulations.csv'), row.names=FALSE)
p <- ggplot(sim, aes(x=generation,y=proportion,group=replicate)) +
  geom_line(alpha=.24,linewidth=.32,color='#317A96') +
  facet_grid(gamma ~ n, labeller=label_both) +
  coord_cartesian(ylim=c(0,1)) + theme_bw(base_size=11) +
  labs(title='Wright-Fisher: 50 replicates per parameter combination',
       x='Generation',y='Cancer cell proportion')
ggsave(file.path(out,'02_WF_trajectories.pdf'),p,width=13,height=9)
wf_stats <- aggregate(proportion ~ n+gamma+generation, data=sim,
                      FUN=function(z) c(mean=mean(z),sd=sd(z)))
wf_stats <- data.frame(n=wf_stats$n, gamma=wf_stats$gamma,
                       generation=wf_stats$generation,
                       mean=wf_stats$proportion[,'mean'],
                       sd=wf_stats$proportion[,'sd'])
write.csv(wf_stats,file.path(out,'03_WF_summary.csv'),row.names=FALSE)

