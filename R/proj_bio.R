#' Projected biomass for 2 years
#'
#' @param report RTMB report object
#' @param Tproj number of years to project, default: 2
#' @export
proj_bio <- function(report, Tproj=2) {

  # values
  F40 = report$F40
  F35 = report$F35
  B40 = report$B40
  Nat = report$Nat
  Sat = report$Sat
  ages = report$ages
  years = report$years
  waa = report$waa
  wt_mature = report$wt_mature
  spawn_frac = report$spawn_fract
  yield_ratio = report$yield_ratio
  M = report$M
  pred_rec = report$pred_rec
  stdev_rec = report$stdev_rec
  A = nrow(Nat) # number of ages
  T = ncol(Nat) # number of years
  slx = report$slx_fish[,T]

  # storage
  N = Cat = Cat_ofl= Zabc = Zofl = S = matrix(0, A, Tproj)
  tot_bio = spawn_bio = F40_proj = F35_proj= rep(0, Tproj)
  # setup
  F40_proj[1] = F40
  F35_proj[1] = F35

  # total F
  Fabc_tot = slx * F40_proj[1]
  Fofl_tot = slx * F35_proj[1]

  # first projection year
  N[1,] = pred_rec
  for(a in 1:(A-1)) {
    N[a+1,1] = Nat[a,T] * Sat[a,T]
  }
  N[A,1] = Nat[A-1,T] * Sat[A-1,T] + Nat[A,T] * Sat[A,T]
  spawn_bio[1] = sum(N[,1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)

  for(t in 1:Tproj) {
    # tier check
    if((spawn_bio[t] / B40) > 1) {
      F40_proj[t] = F40
      F35_proj[t] = F35
    } else {
      F40_proj[t] = F40_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
      F35_proj[t] = F35_proj[t] * (spawn_bio[t] / B40 - 0.05) / 0.95
    }
    # update Fs
    Fabc_tot = slx * F40_proj[t]
    Fofl_tot = slx * F35_proj[t]
    Z = Fabc_tot + M
    Zofl = Fofl_tot + M
    S = exp(-Z)

    # catch
    Cat[,t] = yield_ratio * N[,t] * Fabc_tot / Z * (1 - S)
    Cat_ofl[,t] = yield_ratio * N[,t] * Fofl_tot / Zofl * (1 - exp(-Zofl))

    if(t<Tproj) {
      for(a in 1:(A-1)){
        N[a+1,t+1] = N[a,t] * exp(-yield_ratio * Fabc_tot[a] - M)
      }
      N[A,t+1] = N[A-1,t] * exp(-yield_ratio * Fabc_tot[A-1] - M) +
        N[A,t] * exp(-yield_ratio * Fabc_tot[A] - M)

      tot_bio[t+1] = sum(N[,t+1] * waa)
      spawn_bio[t+1] = sum(N[,t+1] * exp(-yield_ratio * Fabc_tot - M)^spawn_frac * wt_mature)
    }
  }
  catch = colSums(Cat * waa / yield_ratio)
  catch_ofl = colSums(Cat_ofl * waa / yield_ratio)
  tot_bio = colSums(N * waa)

  data.frame(year = max(years)+1:Tproj,
             spawn_bio = spawn_bio,
             tot_bio = tot_bio,
             catch_abc = catch,
             catch_ofl = catch_ofl,
             F40 = F40_proj,
             F35 = F35_proj)
}
