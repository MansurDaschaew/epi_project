library(deSolve)

library(tidyverse)

beta_exp <- function(beta_zero, phi, q, t) {
  beta <- beta_zero * ((1 - phi) * exp(-q * t) + phi)
}

# beta = beta_zero, falls phi = 0 und q = 0.
# vier Variabeln
beta_harm <- function(beta_zero, phi, q, nu, t) {
  beta <- beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)
}

# beta = beta_zero, falls phi = 0, q = 0 und nu = 1.
# vier Variabeln
beta_hyper <- function(beta_zero, phi, q, nu, t) {
  beta <- beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)
}

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
  method <- params["method"]
  N <- S+I+R
  ## now code the model equations
  if(method == 1){
    dSdt <- -beta * S * I/N
    dEdt <- beta * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I
    dRdt <- gamma * I
  }else if(method == 2){
    dSdt <- -beta_exp(beta_zero, phi ,q ,t) * S * I/N
    dEdt <- beta_exp(beta_zero, phi ,q ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I
    dRdt <- gamma * I 
  }else if(method == 3){
    dSdt <- -beta_harm(beta_zero, phi ,q, nu ,t) * S * I/N
    dEdt <- beta_harm(beta_zero, phi ,q, nu ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I
    dRdt <- gamma * I 
  }else if(method == 4){
    dSdt <- -beta_hyper(beta_zero, phi ,q, nu ,t) * S * I/N
    dEdt <- beta_hyper(beta_zero, phi ,q,nu ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I
    dRdt <- gamma * I 
  } else{ stop("Please choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper ")}
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dIdt,dRdt)
  ## return result as a list!
  list(dxdt)
}

#parameters
beta_zero <- 0.8
phi <- 0.4
q <- 0.05
nu <- 0.5

  #beta can be ignored if method is anything else than 1;
  #choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper.

params <- c(alpha=0.5, beta=1, gamma=0.2, method= 4)

times <- seq(from=0,to=60,by=1/60)
xstart <- c(S=9998,E = 1, I=1,R=0)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

#print the output from above (solution of the ode with the given parameters)
out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')