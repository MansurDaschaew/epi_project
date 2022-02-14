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

# In General
## Mittlere Latenzzeit: alpha^-1
## Mittlere infektiöse Zeit: gamma^-1
## Übertragungsrate: beta (R = beta/gamma => beta = gamma * R)
## (phi = gamma/beta_0)


# For Delta
## Alpha = (5,8 - 1,8 - 2)^-1 = 1/2
## Gamma = 1/12
## Beta = 5.5 * 1/12 

beta_zero <-  5.5 / 12 
phi <- 1/(5.5)

# Referenzwerte Xi'an
## erste festgestellte Infektion an Tag x (9.12.21)
## 63 Fälle an an Tag x + 13 (22.12.21)
## Lockdown an Tag x + 14 (23.12.21)

params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.2, method= 1)

times <- 0:30

#times <- seq(from=0,to=60,by=1/60)

xstart <- c(S=12999999,E = 0, I=1,R=0,D=0)


find_delta <- function(deltas,t_end,tol,tol_t=0.05,tol_c=0.1){
  times <- 0:t_end
  results <- data.frame(time = times)
  
  for(i in deltas){
    params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=i, method= 1)
    
    ode(
      func=closed.seir.model,
      y=xstart,
      times=times,
      parms=params
    ) %>%
      as.data.frame() -> out

    results <- cbind(results,out[,'D'])

  }

  head <- c("time",as.character(deltas))
  colnames(results) <- head
  
  help <- as.character(deltas)
  help2 <- rep("time",length(deltas))
  help2 <- paste(help2,help)
  head2 <- rep("",2*length(deltas))

  for(k in 1:length(deltas)){
    head2[2*k - 1] <- help2[k]
    head2[2*k] <- help[k]
  }
  
  dummy <- rep(NA,length.out=length(times))
  selected_results <- dummy

  for(l in 2:(2*length(deltas))){
    selected_results <- cbind(selected_results,dummy)

   }
  colnames(selected_results) <- head2

  for(j in 2:ncol(results)){
    selected <- results[,j]
    relevant <- which(selected>=tol)
    selected <- selected[relevant]
    selected <- results[relevant,c(1,j)]
    l <- nrow(selected[1])
    s1 <- unlist(selected[1],use.names = FALSE)
    s2 <- unlist(selected[2],use.names = FALSE)
    k <- j-1
    k1 <- 2*k - 1
    k2 <- 2*k
    selected_results[1:l,k1] <- s1
    selected_results[1:l,k2] <- s2
  }

  compare <- seq(from=1,by=2,length.out=length(deltas))
  compare <- unlist(selected_results[1,compare],use.names = FALSE)
  compare <- compare[!is.na(compare)]
  compare2 <- compare + 13
  compare3 <- compare + 19
  compare <- c(compare,compare2,compare3)
  c1 <- rep(1,times=length(deltas))
  c2 <- rep(63,times=length(deltas))
  c3 <- rep(175,times=length(deltas))
  comp_time <- c(c1,c2,c3)
  tol_t1 <- tol_t * 63
  tol_t2 <- tol_t * 175
  tol_c1 <- tol_c * 13
  tol_c2 <- tol_c * 19
  comp_time_low <- c(c1,c2- tol_t1,c3-tol_t2)
  comp_time_up <- c(c1,c2+ tol_t1,c3+tol_t2)

  col=rainbow(length(deltas))
  plot(compare,comp_time,ylab="cases",xlab="time",col=col,pch=19,lwd=2,main="Fallzahlen nach Teststrategie und -qualität", sub=paste("Zählgrenze = ",as.character(tol)))
  legend("bottomright",legend=help,pch=19,col=col)
  for(l in 1:length(deltas)){

    lines(rep(compare2[l],times=2),c(c2[l]-tol_t1,c2[l]+tol_t1),col=col[l],lwd=2)
    lines(rep(compare3[l],times=2),c(c3[l]-tol_t2,c3[l]+tol_t2),col=col[l],lwd=2)
  }
  for(l in 1:length(deltas)){
    
    lines(c(compare2[l]-tol_c1,compare2[l]+tol_c1),rep(c2[l],times=2),col=col[l],lwd=2)
    lines(c(compare3[l]-tol_c2,compare3[l]+tol_c2),rep(c3[l],times=2),col=col[l],lwd=2)
  }
  for(l in 1:length(deltas)){
      v1 <- unlist(selected_results[,2*l-1],use.names = FALSE)
      v1 <- v1[!is.na(v1)]
      v2 <- unlist(selected_results[,2*l])
      v2 <- v2[!is.na(v2)]
      lines(v1,v2,col=col[l],lwd=2)
  }

  selected_results
}

a_0.5 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,0.5)
a_0.9 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,0.9)
a_1 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,1)
a_1.1 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,1.1)
a_1.3 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,1.3)
a_1.5 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,1.5)
a_2 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,2)

## Wahl: delta=0.01
a_2.5 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,2.5)

a_3 <- find_delta(c(0.01,0.02,0.03,0.04,0.05),50,3)

# b_0.5 <- find_delta(c(0.06,0.07,0.1,0.5),50,0.5)
# b_0.9 <- find_delta(c(0.06,0.07,0.1,0.5),50,0.9)
# b_1 <- find_delta(c(0.06,0.07,0.1,0.5),50,1)
# b_1.1 <- find_delta(c(0.06,0.07,0.1,0.5),50,1.1)
# b_1.3 <- find_delta(c(0.06,0.07,0.1,0.5),50,1.3)

#####################################################################
#
# a) Kein Eingrif
#
#####################################################################

times <- 0:36
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 1)
xstart <- c(S=12999999,E = 0, I=1,R=0,D=0)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')


find <- out[,'D']
find <- find[which(find >=1)]
found <- sum(find)
found

max(out[,'I'])
max(out[,'E'])
max(out[,'D'])
max(out[,'R'])
min(out[,'S'])/13000000

out[which.min(out$S),]
out[which.max(out$I),]
out[which.max(out$E),]
out[which.max(out$D),]
out[which.max(out$R),]

out[160:200,]
out[200:300,]
out[300:400,]
out[400:500,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]
selection[36,]
selection[37,]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = phi * beta 
# 
# beta = phi * beta_zero = 1/5.5 * 5.5/12 = 1/12
#
#####################################################################

times <- 0:900
params <- c(alpha=0.5, beta=(1/12), gamma=(1/12), delta=0.01, method= 1)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

a %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')
typeof(out)
#outta <- unlist(out[,'I'],use.names = FALSE)

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[1:100,]
out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]
selection[36,]
selection[37,]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_exp
#
# beta_zero * ((1 - phi) * exp(-q * t) + phi
#
# q = 0.05
# 
#####################################################################

q <- 0.05

times <- 0:800
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 2)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)


#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')
#typeof(out)
#outta <- unlist(out[,'I'],use.names = FALSE)
out

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[160:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]
selection[36,]
selection[37,]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_exp
#
# beta_zero * ((1 - phi) * exp(-q * t) + phi
#
# q = 0.005
# 
#####################################################################

q <- 0.005

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 2)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)


#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')


#####################################################################
#
# Variiere beta: beta = beta_exp
#
# beta_zero * ((1 - phi) * exp(-q * t) + phi
#
# q = 0.01
# 
#####################################################################

q <- 0.01

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 2)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)


#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_exp
#
# beta_zero * ((1 - phi) * exp(-q * t) + phi
#
# q = 0.02
# 
#####################################################################

q <- 0.02

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 2)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)


#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_exp
#
# beta_zero * ((1 - phi) * exp(-q * t) + phi
#
# q = 0.05
# 
#####################################################################

q <- 0.05

times <- 0:900
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 2)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)


#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')



#####################################################################
#
# Variiere beta: beta = beta_harm
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)
#
# nu = 1
#
# q = 0.005
#
#####################################################################

nu <- 1
q <- 0.005

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 3)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_harm
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)
#
# nu = 1
#
# q = 0.01
#
#####################################################################

nu <- 1
q <- 0.01

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 3)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_harm
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)
#
# nu = 1
#
# q = 0.02
#
#####################################################################

nu <- 1
q <- 0.02

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 3)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_harm
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)
#
# nu = 1
#
# q = 0.05
#
#####################################################################

nu <- 1
q <- 0.05

times <- 0:400
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 3)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')



#####################################################################
#
# Variiere beta: beta = beta_hyper
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)
#
# nu = 0.5
#
# q = 0.005
#
#####################################################################

nu <- 0.5
q <- 0.005

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 4)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_hyper
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)
#
# nu = 0.5
#
# q = 0.01
#
#####################################################################

nu <- 0.5
q <- 0.01

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 4)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_hyper
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)
#
# nu = 0.5
#
# q = 0.02
#
#####################################################################


#####################################################################
#####################################################################
###
### Hier kommen negative Wert für E und I ab 136 raus...
###
#####################################################################
#####################################################################

nu <- 0.5
q <- 0.02

times <- 0:300
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 4)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Variiere beta: beta = beta_hyper
#
# beta = beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)
#
# nu = 0.5
#
# q = 0.05
#
#####################################################################


#####################################################################
#####################################################################
###
### Hier kommen negative Wert für E und I ab 359 raus...
###
#####################################################################
#####################################################################

nu <- 0.5
q <- 0.05

times <- 0:500
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 4)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[100:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]

#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

#####################################################################
#
# Verschärftes Testkonzept
#
#####################################################################

times <- 0:500
params <- c(alpha=0.5, beta=(5.5/12), gamma=(1/12), delta=0.01, method= 1)
xstart <- c(S=12996450,E = 1097, I=1731,R=666,D=56)

#solving the ode
ode(
  func=closed.seir.model,
  y=xstart,
  times=times,
  parms=params
) %>%
  as.data.frame() -> out

out %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')

out[which.max(out$I),c('time','I')]
out[which.max(out$E),c('time','E')]
out[which.max(out$D),c('time','D')]
out[which.max(out$R),c('time','R')]

out[160:200,]
out[200:300,]
out[300:400,]
out[400:500,]
out[500:600,]
out[600:700,]
out[700:800,]
out[800:900,]
out[900:1000,]
out[1000:1100,]
out[1100:1200,]


#print the output from above (solution of the ode with the given parameters)
selection <- out[,c('time','E','I','D','R')]

selection %>%
  gather(variable,value,-time) %>%
  ggplot(aes(x=time,y=value,color=variable))+
  geom_line(size=2)+
  theme_classic()+
  labs(x='time (days)',y='number of individuals')


