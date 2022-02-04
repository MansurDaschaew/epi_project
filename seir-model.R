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
  D <- x[5]
  ## now extract the parameters
  alpha <- params["alpha"]
  beta <- params["beta"]
  gamma <- params["gamma"]
  delta <- params ["delta"]
  method <- params["method"]
  N <- S+I+R+D+E
  ## now code the model equations
  if(method == 1){
    dSdt <- -beta * S * I/N
    dEdt <- beta * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I - delta * I
    dRdt <- gamma * I + gamma * D
    dDdt <- delta * I - gamma * D
  }else if(method == 2){
    dSdt <- -beta_exp(beta_zero, phi ,q ,t) * S * I/N
    dEdt <- beta_exp(beta_zero, phi ,q ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I - delta * I
    dRdt <- gamma * I + gamma * D
    dDdt <- delta * I - gamma * D
  }else if(method == 3){
    dSdt <- -beta_harm(beta_zero, phi ,q, nu ,t) * S * I/N
    dEdt <- beta_harm(beta_zero, phi ,q, nu ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I - delta * I
    dRdt <- gamma * I + gamma * D
    dDdt <- delta * I - gamma * D
  }else if(method == 4){
    dSdt <- -beta_hyper(beta_zero, phi ,q, nu ,t) * S * I/N
    dEdt <- beta_hyper(beta_zero, phi ,q,nu ,t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I - gamma * D
    dRdt <- gamma * I + gamma * D
    dDdt <- delta * I - gamma * D
  } else{ stop("Please choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper ")}
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dIdt,dRdt,dDdt)
  ## return result as a list!
  list(dxdt)
}

#parameters
beta_zero <- 0.4
phi <- 0.4
q <- 0.05
nu <- 0.5

  #beta can be ignored if method is anything else than 1;
  #choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper.
# Some tests
## 1
params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.4, method= 1)
times <- 0:800

## 2
params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.1, method= 1)
times <-0:400

params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.1, method= 2)
times <- 0:600

#times <- seq(from=0,to=60,by=1/60)


xstart <- c(S=999999,E = 0, I=1,R=0,D=0)

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

