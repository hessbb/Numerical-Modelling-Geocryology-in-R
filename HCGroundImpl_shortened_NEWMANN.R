#== Settings ==================
dt <- 3600    #time step [s]
ns <-  40  # how many time steps to compute before saving?
st <- 40 #365 # length of run [days]
# dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]
zmax  <- 5  # depth of lowermost node center ("how large is the domain") [m]
base  <- 1  # resolution reduction ("how strong is the grid corsened with depth")
# 1: equal spacing, >1 growing
dzmin <- 1
Ti <- -3
# initial temperature [C]

type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity
unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization
bcutype       <- "DIRICHLET"# upper boundary condition type
bcltype       <- "NEWMANN"  # lower boundary condition type
layers.sno    <- c(0)       # index of nodes that are snow
type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity

#== Preparation ================
st    <- st * 3600 * 24 # convert to [s]
t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

#make boundary condition
#bc.upper <- rep(-1.0, nt)
bc.upper <- seq(from = -1, to = 3, length.out = nt) #linear warming
bc.lower <- rep(0, nt)

#make material properties and discretization
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # with of node [m]
                   Tj = rep( Ti, nz), # initial temperature [C]
                   soi = rep(0, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(0, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(2e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

#== Target arrays to hold results =========
out.Tj  <- rbind(NULL,t(mat$Tj))
out.liq <- rbind(NULL,t(mat$liq))
out.t   <- t[0]
out.dz  <- mat$dz

#== Loop over time and solve =======
for (i in 1:nt) {
  mat <- HcGroundImplicit(dt, mat, bc.lower[i], bc.upper[i], unfrozen.type = unfrozen.type,
                          unfrozen.par = unfrozen.par, bcutype = bcutype, bcltype = bcltype,
                          layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
  #output only in specified interval
  check <- i/ns - floor(i/ns)
  if (check == 0) {
    out.Tj  <- rbind( out.Tj,t(mat$Tj))
    out.liq <- rbind(out.liq,t(mat$liq))
    out.t   <- t[i]
    print(paste(i,nt))
  }
}

old.par <- par(mfrow=c(1, 2))

plot(out.Tj[1,],mat$zj, main="Temperature",type="l",xlim=c(-6 , 5))
#== Plotting ================
for (n in 2:147) {
  
  lines(out.Tj[n,], mat$zj, lty = 1, col = "red")
  abline(v=0)
  
}

dev.new("liquid water content")

plot(out.liq[1,],mat$zj, main="Liquid COntent",type="l",xlim=c(0,100))
#== Plotting ================
for (n in 2:147) {
  
  lines(out.liq[n,], mat$zj, lty = 1, col = "red")
  abline(v=0)
  
}

par(old.par)