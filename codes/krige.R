##################################################################################
# R code to read station k0 data from Van Houtte et al. (2018), "A continuous 
# near-surface S-wave attenuation map of New Zealand", Geophysical Journal 
# International (submitted), and calculates the continuous k0 maps using kriging. 
# Please refer to the article, once published, for explanations of 'Model 1' and
# 'Model 2'.
#
# Requires R packages {geoR} and {pracma}.
#
##################################################################################
# Read data:
# - Coordinates are NZTM northings and eastings, converted to km for convenience.
# - Data are log10(k0)
k0_dat <- readRDS("nz_k0_data.rds")
# Initial values for regression
ini.nug <- 0
ini.phi <- 100
ini.sill <- 0.1
# Matern order (fixed)
theta <- 0.5 
# Model 1 - solving for nugget
model1 <- likfit(k0_dat, cov.model="matern", ini.cov.pars = c(ini.sill, ini.phi), 
                 nug = ini.nug, fix.nugget=F, kappa=theta, fix.kappa=T,trend=~TVZ)
# Model 2 - fixing the nugget
fixed.nugget <- log10(1.5)^2
model2 <- likfit(k0_dat, cov.model="matern", ini.cov.pars = c(ini.sill,ini.phi), 
                 nug = fixed.nugget, fix.nugget=T, kappa=theta, fix.kappa=T, 
                 trend=~TVZ)

# Prediction grid in NZTM coordinates, i.e. northing and easting, converted to km
pred.grid<- expand.grid(seq(2000000, 3000000, l=1000), 
                        seq(5200000, 6700000, l=1000)) / 1000
# Read digitised 'whole TVZ' model of Wilson et al. (1995), coordinates in km
whole.tvz=read.table("whole_tvz_ne.txt", header = TRUE)
ftvz<-ifelse(inpolygon(pred.grid$Var1, pred.grid$Var2, 
                       whole.tvz$Easting/1000, whole.tvz$Northing/1000), 1, 0)
# Prediction
nug.model1 <- model1$tausq
theta.model1 <- 0.5
sill.model1 <- model1$sigmasq
phi.model1 <- model1$phi
psiA.model1 <- model1$aniso.pars[1]
psiR.model1 <- model1$aniso.pars[2]

kc.model1 <- krige.conv(k0_dat, loc = pred.grid, 
                        krige = krige.control(type.krige="ok", cov.model="matern",
                                              cov.pars=c(sill.model1, phi.model1), 
                                              kappa=theta.model1, nugget=nug.model1,
                                              aniso.pars=c(psiA.model1, psiR.model1), 
                                              trend.d=~TVZ, trend.l=~ftvz),
                        output = output.control(signal=T))

nug.model2 <- model2$tausq
theta.model2 <- 0.5
sill.model2 <- model2$sigmasq
phi.model2 <- model2$phi
psiA.model2 <- model2$aniso.pars[1]
psiR.model2 <- model2$aniso.pars[2]

kc.model2 <- krige.conv(k0_dat, loc = pred.grid, 
                        krige = krige.control(type.krige="ok",cov.model="matern",
                                              cov.pars=c(sill.model2, phi.model2), 
                                              kappa=theta, nugget=nug.model2,
                                              aniso.pars=c(psiA.model2,psiR.model2),
                                              trend.d=~TVZ,trend.l=~ftvz), 
                        output = output.control(signal=T))

# Units of grid back to m
Nor <- kc.model1$Var2 * 1000
Eas <- kc.model2$Var1 * 1000
k0.model1 <- 10^kc.model1$predict
k0.model2 <- 10^kc.model2$predict
