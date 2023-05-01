library("data.table")
library("ggplot2")
library("blockrand")
library("rstan")
library("poisson")

# Same parameterisation as used in stan model
emax2_p <- function(dose, e0, emax_delta, ed50, hill, prob = T){
  tau <- hill * ( log(ed50) - log(dose[-1]) )
  nu <- -log(1 + exp(tau))
  if(hill < 0){ eta <- e0 + emax_delta * c(1, exp(nu)) }
  if(hill > 0){ eta <- e0 + emax_delta * c(0, exp(nu)) }
  if(hill == 0){ stop("zero hill not supported") }
  if(prob == T){ plogis(eta) } else { eta }
}

deriv_emax2_p <- function(dose, e0, emax_delta, ed50, hill){
  # prevent infinite (but not using the zero dose anyway)
  dose[dose == 0] <- .Machine$double.eps^12
  u <- e0 * dose^hill + e0 * ed50^hill + emax_delta * dose^hill
  v <- dose^hill + ed50^hill
  up <- hill * dose^(hill-1) * (e0 + emax_delta)
  vp <- hill * dose^(hill-1)
  der = (up * v - u * vp)/(v^2)
  der
}

get_study_cohort <- function(
    lpar
){
  budget_remaining <- lpar$budget_tot * ((100-lpar$budget_buffer_pct)/100)
  i <- 1
  ld <- list()
  while(budget_remaining > 0){

    # create a block
    d_alloc <- data.table(
      blockrand(
        1, num.levels = length(lpar$dose), levels = seq_along(lpar$dose),
        # resulting block sizes is 2 * 1:5 = {2, 4, 6, 8, 10}
        block.sizes = 1:2,
        block.prefix = paste0("B")
      ))
    names(d_alloc) <- c("id", "id_block", "blocksize", "id_arm")

    d_alloc[, nom_gp := rbinom(.N, 1, lpar$prob_gp_nom)]
    d_alloc[, dose := lpar$dose[id_arm]]
    d_alloc[, p_k := lpar$p_k[id_arm]]
    d_alloc[, y := rbinom(.N, 1, p_k)]
    # only those that initiate trt go on to svr
    d_alloc[y == 1, svr := rbinom(.N, 1, lpar$prob_pt_svr)]
    d_alloc[y == 0, svr := 0]

    # Costings just proxy, might need reworking.
    d_alloc[y==1, cost := lpar$pt_admin + # all recv admin
              1000 * dose +               # success recv dose
              nom_gp * lpar$gp_coincen +  # success recv dose
              lpar$gp_admin +             # all recv admin
              lpar$pt_admin_svr * svr     # success recv svr admin if fu
    ]
    d_alloc[y==0, cost := lpar$pt_admin + # all recv admin
              lpar$gp_admin               # all recv admin
    ]
    ld[[i]] <- copy(d_alloc)
    i <- i + 1
    budget_remaining <- budget_remaining - sum(d_alloc$cost)
  }
  d <- rbindlist(ld)
  d[, `:=`(id_block=NULL,blocksize=NULL)]
  tt <- c(0, nhpp.event.times(lpar$nhpoislambda, nrow(d)-1, lpar$nhomogintens))
  d[, t0 := tt]
  d[, id := 1:.N]
  d
}

main_util <- function(){

  lpar = list(
    # Dose defined globally at head of file - factor out shortly.
    e0 = qlogis(0.2),
    emax_delta = qlogis(0.3) - qlogis(0.2),
    ed50 = 400,
    hill = 3,
    dose = seq(0, 1, len = 21), # i.e. 0 to 1000 by 50
    p_k = NULL, # emax2_p(dose, e0, emax_delta, ed50, hill),
    med_delta = 0.1,
    max_imbal = 5,
    # total budget available
    budget_tot = 700000,
    # interim size is based on budget expenditure
    budget_intrm = 50000,
    # but interim 1 has a minimum fixed size
    n_interim_1 = 500,
    budget_buffer_pct = 10,
    # percentage that nominate gp
    prob_gp_nom = 0.5,
    prob_pt_svr = 0.5,
    pt_admin = 105,
    pt_admin_svr = 75,
    gp_admin = 50,
    gp_coincen = 100,
    # failure time 84 days = 12 wks
    t_fail = 84,
    # min time to success 14 days
    t_min_success = 14,
    # enrolment times - poisson process
    # events per day
    nhpoislambda = 780/365,
    # ramp up over 6 months (180 days)
    nhomogintens = function(t) pmin(t/180, 1),
    debug = T,
    debug_err = F,
    mcmciter = 7000
  )
  # refine as necessary
  lpar$e0 <- qlogis(0.15)
  lpar$emax_delta <- qlogis(0.45) - qlogis(0.15)
  lpar$ed50 <- 0.3
  lpar$hill <- 2

  lpar$p_k <- emax2_p(
    lpar$dose, lpar$e0, lpar$emax_delta, lpar$ed50, lpar$hill)

  K <- length(lpar$p_k)


  d_all <- get_study_cohort(lpar)
  d_mod <- d_all[, .(y = sum(y), n = .N), keyby = dose]
  # To preserve indexing in the models, we want the full design space
  # (all possible doses) to be represented in the data even if they got
  # no allocation nor had any events.
  d_mod <- merge(data.table(dose = lpar$dose),d_mod,by = "dose",all.x = T)
  d_mod[is.na(y), `:=`(y =0, n = 0)]

  # quick plot truth versus observed
  plot(lpar$dose, lpar$p_k, type = "l")
  points(lpar$dose, d_mod$y/d_mod$n)

  lsdat = list(
    N = nrow(d_mod), y = d_mod[, y], n = d_mod[, n], dose = d_mod$dose,
    prior_only = F, debug = F, pri_mu = c(0, 0), pri_s = c(2.5,1,1,0.85),
    pri_nu = 4
  )

  # Model
  mod1 <- rstan::stan_model("./stan/emax2.stan")
  f1 <- rstan::sampling(
    mod1,
    data = lsdat, iter = lpar$mcmciter, refresh = 0, chains = 1
  )

  # posterior
  post_emax <- as.matrix(f1, pars = c("b0","bemax","bed50","bhill"))
  colMeans(post_emax)
  post_p <- as.matrix(f1, pars = paste0("p"))
  colMeans(post_p)

  # plot probs
  dfig <- melt(data.table(post_p), measure.vars  = colnames(post_p))
  dfig[, id_arm := gsub("p[", "", variable, fixed = T)]
  dfig[, id_arm := as.integer(gsub("]", "", id_arm, fixed = T))]
  dfig <- dfig[, .(p_mu = mean(value),
           q_025 = quantile(value,0.025),
           q_975 = quantile(value,0.975)), keyby = id_arm]
  setkey(dfig, id_arm)
  dfig[, p_tru := lpar$p_k]
  ggplot(dfig, aes(x = id_arm, y = p_mu)) +
    geom_linerange(aes(ymin = q_025, ymax = q_975)) +
    geom_point() +
    geom_line(aes(x = id_arm, y = p_tru))

  # plot emax
  dfig <- melt(data.table(post_emax), measure.vars  = colnames(post_emax))
  dfig <- dfig[, .(mu = mean(value),
                   q_025 = quantile(value,0.025),
                   q_975 = quantile(value,0.975)), keyby = variable]
  dfig[variable == "b0", tru := lpar$e0]
  dfig[variable == "bemax", tru := lpar$emax_delta]
  dfig[variable == "bed50", tru := lpar$ed50]
  dfig[variable == "bhill", tru := lpar$hill]
  ggplot(dfig, aes(x = variable, y = mu)) +
    geom_linerange(aes(ymin = q_025, ymax = q_975)) +
    geom_point() +
    geom_point(aes(x = variable, y = tru), col = 2)

  # quantities of interest
  post_p_minus_p0 <- sapply(2:ncol(post_p), function(ii){
    post_p[,ii] - (post_p[, 1] + lpar$med_delta)
  })
  pr_abs_dist_med <- prop.table(
    table(factor(max.col(-abs(post_p_minus_p0))+1,levels = 2:K,labels = 2:K)))
  pr_ed <- apply(post_p_minus_p0, 2, function(z){ mean(z>0) })
  eta_der <- do.call(rbind, lapply(1:nrow(post_emax), function(ii){
    deriv_emax2_p(lpar$dose,
                  post_emax[ii,1],
                  post_emax[ii,2],
                  post_emax[ii,3],
                  post_emax[ii,4])
  }))
  colMeans(eta_der)[-1]

}

