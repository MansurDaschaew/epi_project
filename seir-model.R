library(deSolve)

library(tidyverse)



closed.seir.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  E <- x[2]
  I <- x[3]
  R <- x[4]
  ## now extract the parameters
  alpha <- params["alpha"]
  beta <- params["beta"]
  gamma <- params["gamma"]
  N <- S+I+R
  ## now code the model equations
  dSdt <- -beta * S * I/N
  dEdt <- beta * S * I/N - alpha * E 
  dIdt <- alpha * E - gamma * I
  dRdt <- gamma * I
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

parms <- c(alpha=0.5, beta=0.5, gamma=0.2)

times <- seq(from=0,to=60,by=1/60)
xstart <- c(S=9998,E = 1, I=1,R=0)

ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=parms
) %>%
  as.data.frame() -> out
  
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')