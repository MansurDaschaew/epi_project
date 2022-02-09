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

#oszilierend beta 
#1 konstant oszillierend z.b. jahreszeiten abhängig
beta_cos <- function(beta_zero, t, c1 = 0.25, c2 = 365) {
  beta <- beta_zero + cos(t*2*pi/c2) * c1
  return(beta)
}
#2 ansteigend z.b jahreszeiten abhängig mit ansteigender pandemiemüdigkeit 
#3 absteigend (notwendigkeit??)
#ansteigende beta 
expdata <- data.frame(a0 = 0:100, a1 = beta_exp(0.5, 1.2, 0.1,0:100))
ggr <- ggplot() + geom_line(data = expdata, aes(a0, a1),size =1)



expadata <- data.frame(a0 = 1:365, a1 = expy(beta_zero = 1, x = 1:365))
ggra <- ggplot() + geom_line(data = expadata, aes(a0, a1),size =1)

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
  } else if(method == 5){
    dSdt <- -expy(beta_zero = beta_zero, x = t) * S * I/N
    dEdt <- expy(beta_zero = beta_zero, x = t) * S * I/N - alpha * E 
    dIdt <- alpha * E - gamma * I - gamma * D
    dRdt <- gamma * I + gamma * D
    dDdt <- delta * I - gamma * D
  }else{ stop("Please choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper ")}
  ## combine results into a single vector
  dxdt <- c(dSdt,dEdt,dIdt,dRdt,dDdt)
  ## return result as a list!
  list(dxdt)
}

#parameters
beta_zero <- 0.3
phi <- 0.8
q <- 0.05
nu <- 0.5

  #beta can be ignored if method is anything else than 1;
  #choose method= 1 for const beta, method= 2 for beta_exp, method= 3 for beta_harm, method= 4 for beta_hyper.
# Some tests
## 1
params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.4, method= 1)
times <- 0:800

params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.1, method= 2)
times <- 0:600

#times <- seq(from=0,to=60,by=1/60)

## 2
expy <- function(beta_zero,c1=0.06,c2=0.009,x){
  beta <- beta_zero*((c1*abs(cos(c2*pi*x))*(log(x)))+1)
  return(beta)
}
expy <- function(beta_zero,c1=0.05,c2=0.01,x){
  beta <- beta_zero*((c1*(cos(c2*pi*x))^2*(log(x)))+1)
  return(beta)
}
params <- c(alpha=0.1, beta=0.7, gamma=0.1, delta=0.1, method= 5)
times <-1:300



xstart <- c(S=99999,E = 0, I=1,R=0,D=0)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

#print the output from above (solution of the ode with the given parameters)
out[c(1,3,4,6)] %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

