library(tibble)

Fin <- 120 #max. simulation time
N <- 200 #population size
gamma <- 0.15 #recovery rate
beta <- 0.2 / (N + 1) #contact rate
lam <- N * (N + 1) / 2 * beta
#results <- []

# beta = beta_zero, falls phi = 0 und q = 0.
# drei Variabeln
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


# 0 < phi < 1? Sonst beta negativ.

# Funktion, die beta_exp, beta_harm und beta_hyper von 0 bis tFin plottet
# Mögliche Argumente für fun_beta sind "exp", "harm" und "hyper"
# Höchstens ein Elemente der Menge {beta_zero, gamma, q, nu} darf mehrdimensional sein.
# phi = gamma / beta_zero

plot_beta <-
  function(tFin,
           fun_beta = "exp",
           beta_zero = 0.1,
           gamma = 0,
           q = 0,
           nu = 1) {
    t <- seq(0, tFin, length.out = 5 * tFin)
    phi <- gamma / beta_zero
    
    if(length(phi[which(phi > 1) >= 1])){
      stop("phi = gamma / beta_zero muss kleiner als 1 sein!")
    }
    
    # Überprüfe, ob mehr als ein Elemente der Menge {beta_zero, gamma, q, nu} mehrdimensional ist.
    # Brich ab, falls das der Fall ist.
    x <- c(length(beta_zero),length(gamma),length(q),length(nu))
    if(length(which(x > 1)) > 1){
      stop("Höchstens ein Argument der Menge {beta_zero, gamma, q, nu} darf mehrdimensional sein!")
    }
    
    # Falls genau ein Elemente der Menge {beta_zero, phi, q, nu} mehrdimentional ist,
    # finde dieses, initialisiere eine Matrix, in die die Werte von beta geschrieben werden.
    # Bestimme anderenfalls nur die Werte zu beta.
    # Bestimme außerdem den Vektor und den Titel für die Legende des Plots.
    if (length(beta_zero) > 1) {
      legend <- beta_zero
      title <- "beta_zero"
      legend3 <- phi
      title3 <- "phi"
      
      beta <-
        matrix(rep(0, times = length(beta_zero) * length(t)), ncol = length(beta_zero))
      
      for (i in 1:length(beta_zero)) {
        if (fun_beta == "exp") {
          beta[, i] <- beta_exp(beta_zero[i], phi[i], q, t)
          
        } else if (fun_beta == "harm") {
          beta[, i] <- beta_harm(beta_zero[i], phi[i], q, nu, t)
          
        } else if (fun_beta == "hyper") {
          beta[, i] <- beta_hyper(beta_zero[i], phi[i], q, nu, t)
        }
      }
    } else if (length(gamma) > 1) {
      legend <- phi
      title <- "phi"
      
      beta <- matrix(rep(0, times = length(phi) * length(t)), ncol = length(phi))
      
      for (i in 1:length(phi)) {
        if (fun_beta == "exp") {
          beta[, i] <- beta_exp(beta_zero, phi[i], q, t)
          
        } else if (fun_beta == "harm") {
          beta[, i] <- beta_harm(beta_zero, phi[i], q, nu, t)
          
        } else if (fun_beta == "hyper") {
          beta[, i] <- beta_hyper(beta_zero, phi[i], q, nu, t)
        }
      }
    } else if (length(q) > 1) {
      legend <- q
      title <- "q"
      
      beta <- matrix(rep(0, times = length(q) * length(t)), ncol = length(q))
      
      for (i in 1:length(q)) {
        if (fun_beta == "exp") {
          beta[, i] <- beta_exp(beta_zero, phi, q[i], t)
          
        } else if (fun_beta == "harm") {
          beta[, i] <- beta_harm(beta_zero, phi, q[i], nu, t)
          
        } else if (fun_beta == "hyper") {
          beta[, i] <- beta_hyper(beta_zero, phi, q[i], nu, t)
        }
      }
    } else if (length(nu) > 1) {
      legend <- nu
      title <- "nu"
      
      beta <- matrix(rep(0, times = length(nu) * length(t)), ncol = length(nu))
      
      for (i in 1:length(nu)) {
        if (fun_beta == "harm") {
          beta[, i] <- beta_harm(beta_zero, phi, q, nu[i], t)
          
        } else if (fun_beta == "hyper") {
          beta[, i] <- beta_hyper(beta_zero, phi, q, nu[i], t)
        }
      }
    } else{
      if (fun_beta == "exp") {
        beta <- beta_exp(beta_zero, phi, q, t)
        
      } else if (fun_beta == "harm") {
        beta <- beta_harm(beta_zero, phi, q, nu, t)
        
      } else if (fun_beta == "hyper") {
        beta <- beta_hyper(beta_zero, phi, q, nu, t)
      }
    }
    
    # Für die zweite Legende (eindimensionale Parameter)
    # Wandele Paramater in character-Vektoren um
    char_q <- paste("q =", as.character(q))
    char_phi <- paste("phi =", as.character(phi))
    char_beta_zero <- paste("beta_zero =", as.character(beta_zero))
    char_nu <- paste("nu =", as.character(nu))
    legend2 <- ""
    
    # Wähle den Text unter dem Plot (die Funktionsvorschrift) und die eindimensionalen Parameter
    if(fun_beta == "exp"){
      sub <- "beta_exp = beta_zero * ((1 - phi) * exp(-q * t) + phi)"
      if(length(beta_zero) == 1){
        legend2 <- c(legend2, char_beta_zero)
      }
      if(length(phi) == 1){
        legend2 <- c(legend2, char_phi)
      }
      if(length(q) == 1){
        legend2 <- c(legend2, char_q)
      }
    }
    if(fun_beta == "harm"){
      sub <- "beta_harm = beta_zero * ((1 - phi) / (1 + q * nu * t) + phi)"
      if(length(beta_zero) == 1){
        legend2 <- c(legend2, char_beta_zero)
      }
      if(length(phi) == 1){
        legend2 <- c(legend2, char_phi)
      }
      if(length(q) == 1){
        legend2 <- c(legend2, char_q)
      }
      if(length(nu) == 1){
        legend2 <- c(legend2, char_nu)
      }
    }
    if(fun_beta == "hyper"){
      sub <- "beta_hyper = beta_zero * ((1 - phi) / (1 + q * nu * t) ^ (1/nu) + phi)"
      if(length(beta_zero) == 1){
        legend2 <- c(legend2, char_beta_zero)
      }
      if(length(phi) == 1){
        legend2 <- c(legend2, char_phi)
      }
      if(length(q) == 1){
        legend2 <- c(legend2, char_q)
      }
      if(length(nu) == 1){
        legend2 <- c(legend2, char_nu)
      }
    }
    legend2 <- legend2[legend2 != ""]
    
    # Erstelle den Plot.
    # Die Vorgehensweise unterscheidet sich je nachdem, ob ein Wert 
    # der Menge {beta_zero, phi, q, nu}  mehrdimensional ist.
    if(is.matrix(beta)) {
      M <- max(beta)
      m <- min(beta)
      col=rainbow(ncol(beta))
      plot(t, beta[,1], type = "l",lwd=2,col=col[1],ylim=c(m,M),ylab="beta",main="Beta", sub=sub)
      legend("topright",legend=legend,col=col,pch=19,title=title)
      legend("bottomleft",legend=legend2)
      if(length(beta_zero) > 1){
        legend("top",legend=legend3,col=col,pch=19,title=title3)
      }
      for(j in 2:ncol(beta)){
        lines(t,beta[,j], lwd=2,col=col[j])
      }
    }else {
      plot(t, beta, type="l")
    }
    
  }

# beta_zero mehrdimensional => phi variiert ebenfalls -> drei Legenden
plot_beta(80, "harm", c(0.3,0.4,0.5), 1/6, 0.01)

# phi >= 1 Fehlermeldung
plot_beta(80, "harm", 0.4, c(0.3,0.4,0.5), 0.01)

# gamma mehrdimensional
plot_beta(80, "harm", 0.4, c(0.15, 0.2, 0.25, 0.3), 0.01)

# Bsp. aus Buch:
# beta_zero = 0.4
# gamma = 1/6
# => phi = gamma / beta_zero = 1/6 * 5/2 = 5/12
# q = (0,0.005,0.01,0.02,0.05)
# harmonic: nu = 1
# hyper: nu = 0.5

plot_beta(80, "exp", 0.4, 1/6, c(0,0.005,0.01,0.02,0.05))

plot_beta(80, "harm", 0.4, 1/6, c(0,0.005,0.01,0.02,0.05),1)

plot_beta(80, "hyper", 0.4, 1/6, c(0,0.005,0.01,0.02,0.05),0.5)

single_run <-
  function(tFin,
           N,
           gamma,
           beta,
           fun_beta = "const",
           beta_zero = 0.1,
           phi = 0,
           q = 0,
           nu = 1) {
    beta <- beta / (N + 1)
    lam <- N * (N + 1) / 2 * beta
    xx <-
      rep(c(0), times = N + 1) #all individuals in state 0 = susceptible
    xx[1] <- 1 #wlog individual 0 as index case
    tt <- c(0)
    t <- tt[1]
    ss <- xx
    St <- which(xx == 0)
    It <- which(xx == 1)
    Rt <- which(xx == 2)
    nI <- length(It)
    nS <- length(St)
    nR <- length(Rt)
    sus <- nS
    inf <- nI
    rec <- nR
    
    stepCount <- 0
    
    while (length(nI) > 0 && t <= tFin) {
      stepCount <- stepCount + 1
      
      if (nS > 0) {
        rate <- nI * gamma + lam
        u <- runif(1)
      } else{
        rate <- nI * gamma
        u <- -1
      }
      pRec <- nI * gamma / rate
      dt <- rexp(1, rate = c(gamma)) #size of time step
      t <- t + dt #update current time
      
      if (u < pRec) {
        ###randomly choose index among infectious individuals
        i <- It[runif(1, 1, nI)]
        ###set status of this individual to 2 = recovered
        xx[i] <- 2
        ###store t as time when something happened
        tt <- c(tt, t)
        St <- which(xx == 0)
        It <- which(xx == 1)
        Rt <- which(xx == 2)
        nI <- length(It)
        nS <- length(St)
        nR <- length(Rt)
        sus <- c(sus, nS)
        inf <- c(inf, nI)
        rec <- c(rec, nR)
      } else{
        ### randomly choose one individual
        i <- as.integer(runif(1, 1, N + 1))
        j <- i
        ### randomly choose another (different) individual
        while (i == j) {
          j <- as.integer(runif(1, 1, N + 1))
        }
        
        if (xx[i] + xx[j] == 1) {
          ### in case one is infectious and the other susceptible
          ### make both infectious
          xx[i] <- 1
          xx[j] <- 1
          ### and store t as time when something happened
          tt <- c(tt, t)
          St <- which(xx == 0)
          It <- which(xx == 1)
          Rt <- which(xx == 2)
          nI <- length(It)
          nS <- length(St)
          nR <- length(Rt)
          sus <- c(sus, nS)
          inf <- c(inf, nI)
          rec <- c(rec, nR)
        }
        
      }
      if (fun_beta == "exp") {
        beta <- c(beta, beta <- beta_exp(beta_zero, phi, q, t))
      } else if (fun_beta == "harm") {
        beta <- c(beta, beta_harm(beta_zero, phi, q, nu, t))
      } else if (fun_beta == "hyper") {
        beta <- c(beta, beta_hyper(beta_zero, phi, q, nu, t))
      }
    }
    
    inf2 <- rep(c(N + 1), times = length(sus)) - sus - rec
    inf_max <- max(inf2)
    sus_end <- sus[length(sus)]
    
    # plot(tt, sus,
    #      type ="l",
    #      col="blue",
    #      main = "Titel des Plots",
    #      ylim = c(1,N),
    #      xlim = c(0,tFin),
    #      xlab = "Time",
    #      ylab = "Individuals")
    # lines(tt, inf,
    #       type="l",
    #       col="red")
    # lines(tt, rec,
    #       type="l",
    #       col="green")
    
    # plot_single <- plot(tt, sus,
    #      type ="l",
    #      col="blue",
    #      main = "Titel des Plots",
    #      ylim = c(1,N),
    #      xlim = c(0,tFin),
    #      xlab = "Time",
    #      ylab = "Individuals")
    # lines(tt, inf2,
    #       type="l",
    #       col="red")
    # lines(tt, rec,
    #       type="l",
    #       col="green")
    
    output <- c(inf_max, sus_end)
  }

multiple_runs <-
  function(tFin,
           N,
           gamma,
           runs,
           beta = "const",
           beta_zero = 0,
           phi = 0,
           q = 0,
           nu = 1) {
    inf_max <-
      rep(c(0), times = runs) # Vektor für die maximale Anzahl an Infizierten pro Durchlauf.
    sus_end <-
      rep(c(0), times = runs) # Vektor für verbleibende empfängliche Individuen pro Durchlauf.
    count_max_inf <-
      rep(c(0), times = runs) # Hilfsvektor zum Zählen der Häufigkeit, dass i=1,...,N
    # die maximale Anzahl an Infizierten ist.
    count_sus_end <-
      rep(c(0), times = runs) # Hilfsvektor zum Zählen der Häufigkeit, dass am Ende
    # i=1,...N, Empängliche verbleiben.
    counter <- c(1:(N + 1))
    counter2 <- c(0:N)
    
    # Führe die Simulationen durch und trage bei jedem Durchlauf die maximale Anzahl an Infizierten und die
    # Anzahl der verbleibenden Empfänglichen in die zuvor initialisierten Vektoren ein.
    for (i in 1:runs) {
      sim <- single_run(tFin, N, gamma, beta)
      
      inf_max[i] <- sim[1]
      sus_end[i] <- sim[2]
    }
    
    # Zähle die maximale Anzahl Infizierter: Im j-ten Eintrag steht, wie oft j die maximale Anzahl an
    # Infizierten war.
    for (j in 1:N + 1) {
      count_max_inf[j] <- length(which(inf_max == j))
    }
    
    # Zähle die Anzahl verbleibender Empfänglicher: Im k-ten Eintrag steht, wie oft k-1 die Anzahl an
    # am Ende übrig gebliebenen Empfänglichen war.
    for (k in 0:N) {
      count_sus_end[k + 1] <- length(which(sus_end == k))
    }
    
    # Wähle die Indizes größer 0 aus. Der Index gibt die Anzahl der verbleibenden Empfänglichen an,
    # der Eintrag die Häufigkeit dieser Anzahl unter allen Simulationen.
    sim_sus <- count_sus_end[which(count_sus_end > 0)]
    # Erhalte die Anzahl der verbleibenden Empfänglichen. Indexkorrektur.
    sus <- (which(count_sus_end > 0) - 1)
    
    # Analog
    sim_inf <- count_max_inf[which(count_max_inf > 0)]
    inf <- which(count_max_inf > 0)
    
    # Berechne statische Größen
    ### Wähle Ausbrüche aus
    
    outbreaks <- N - sus_end[which(sus_end < N)]
    ### Berechne Mittelwerte und Median der Infizierten
    mean_outbreaks <- mean(outbreaks)
    median_outbreaks <- median(outbreaks)
    ### Berechne den Anteil der Ausbrüche unter den Simulationen
    num_outbreaks <- length(which(sus_end < N))
    perc <- num_outbreaks / runs * 100
    
    #Infos auf der Konsole
    print(
      paste(
        "Es wurden",
        runs,
        "Simulationen mit jeweils",
        N,
        "empfänglichen Individuen und",
        tFin,
        "Zeitschritten durchgeführt."
      )
    )
    
    print(
      paste(
        "In",
        num_outbreaks,
        "von",
        runs,
        "Simulationen kam es zum Ausbruch, das entspricht",
        perc,
        "% der Fälle. Bei einem Ausbruch kam es durchschnittlich zu",
        mean_outbreaks,
        "Infektionen. Der Median ist",
        median_outbreaks,
        "."
      )
    )
    
    output <-
      list(
        tibble(inf_max = inf_max, sus_end = sus_end),
        tibble(sus = sus, sim_sus = sim_sus),
        tibble(inf = inf, sim_inf = sim_inf),
        c(mean_outbreaks, median_outbreaks)
      )
  }

plot_sus <- function(x) {
  y <- x[[2]]
  plot(y$sus,
       y$sim_sus,
       type="h",
       lwd="3",
       main = "Vebleibende Empfängliche",
       ylab = "Anzahl der betroffenen Simulationen",
       xlab = "Anzahl der verbleibenden Empfänglichen")
  
}

plot_inf <- function(x) {
  y <- x[[3]]
  z <- x[[4]]
  plot(y$inf,
       y$sim_inf,
       type="h",
       lwd="3",
       main = "Maximale Anzahl an Infizierten",
       xlab = "Maximale Anzahl an Infizierten",
       ylab = "Anzahl der betroffenen Simulationen")
  lines(c(min(y$inf),max(y$inf)),c(z[1],z[1]),col="blue",lty=2)
  lines(c(min(y$inf),max(y$inf)),c(z[2],z[2]),col="green",lty=3)
  legend("topright",legend=c("mean","median"),col=c("blue","green"),pch=19)
}

a <- multiple_runs(50, 100, 0.5, 10000, 0.9, "hyper", 0.9)

plot_sus(a)

plot_inf(a)
