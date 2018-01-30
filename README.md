# Numerical-Modelling-Geocryology-in-R
Here is some code I am working on that will predict how a section of ground will react to different temperature variations over the course of a year, namely looking at ice penetration at depth.
# Analytical-exercise.R
#Original code in collaboration with Robert Martin and using packages created by Stephan Gruber
# Formulas:
#'
#'             CALONNE: Calonne, N., Flin, F., Morin, S., Lesaffre, B., du 
#'             Roscoat, S. R., & Geindreau, C. (2011). Numerical and experimental
#'             investigations of the effective thermal conductivity of snow. 
#'             Geophysical Research Letters, 38(23), doi:10.1029/2011GL049234  
#'             "YEN": Yen, Y.-C. (1981). Review of thermal properties of 
#'              snow, ice and sea ice (34 pages). Hanover, NH, USA.
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             STURM: Sturm, M., J. Holmgren, M. König, and K. Morris (1997),
#'             The thermal conductivity of seasonal snow, Journal of 
#'             Glaciology, 43(143), 26–41.
#'
#'             JORDAN: Jordan, R. E., Andreas, E. L., & Makshtas, A. P. 
#'             (1999). Heat budget of snow-covered sea ice at North Pole 4. 
#'             Journal of Geophysical Research, 104(C4), 7785. 
#'             doi:10.1029/1999JC900011
#'
#'             WILLIAMS: Figure 4.11 in "The Frozen Earth: fundamentals of 
#'             geocryology" by P. Williams and M. Smith, 1989.
#' 
#' @param type.gnd Character string indicating the type of parameterization to be
#'             used for ground; the standard is "COSENZA". Available are:
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             GEOMETRIC: Geometric mean, intermediate mixed conductivity model, 
#'             approximation of randomly oriented consituent elements. 
#'
#'             ARITHMETIC: Arithmetic mean, high mixed conductivity model, 
#'             approximation of consituent elements layered parallel to 
#'             temperature gradient.
#'
#'             HARMONIC: Harmonic mean, low mixed conductivity model, 
#'             approximation of consituent elements layered normal to 
#'             temperature gradient.
#'
#Input Parameters - Make each of these inputs to a function later
Geotherm <- 0.2  
Vol_HC <- 1.8e6 #Volumetric Heat Capacity
Kt <- 2.5 #Thermal Conductivity
Period <- 365*24*3600
omega <- 2*pi/Period #How quickly the amplitude decreases w depth
Amp <- 4 
MAGST <- -0.5
kdiff <- Kt/Vol_HC #How easily the heat can move through the ground thermal diffusivity

n_time_steps=700
n_depth_steps=40

dt=5*3600*24
t=n_time_steps*dt

dz=.5
z=n_depth_steps*dz

delT=rep(0, n_depth_steps*n_time_steps)
dim(delT)=c(n_depth_steps,n_time_steps)

for (z_int in 1:n_depth_steps){
  z=(z_int-1)*dz
  
  Initial_Offset=MAGST + (Geotherm/Kt)*z
  for (t_int in 1:n_time_steps){
    
    t=(t_int-1)*dt
    
    delT[z_int,t_int] <- Initial_Offset +Amp * exp(-z*sqrt(omega/(2*kdiff)))*cos(omega*t-z*sqrt(omega/(2*kdiff)))
  
  }
  
}


plot((1:t_int)/73,delT[1,],type="l",xlab='years',ylim=c(-10,5))

for (depthnumber in 2:n_depth_steps){
  lines((1:t_int)/73,delT[depthnumber,],type="l") 
}
#Numerical1
dt <- 3600*2  #Each Time-step is 2 hours
ns <- 24*5 #Info saved every 10 days
st <- 365 # 10 years required for spin-up (initial experiment was 10 years, wasn't able to determine if Spin-up was done)
zmax  <- 5  # Lowest Depth Modelled
base  <- 1.45  # resolution reduction ("how strong is the grid corsened with depth")
# 1: equal spacing, >1 growing
dzmin <- 0.05 # dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]
Ti <- -15
omega=2*pi/(24*365*3600)
Seasonal_Amp=4
MAGST = -.5

#Difference between MAGST and Ti 14.5
#PartB Question 2
#Experiment 1
dt <- 3600*24*5  #Each Time-step is 5 days
ns <- 1 #Every step is saved
st <- 36500 # 100 years required for spin-up (initial experiment was 10 years, wasn't able to determine if Spin-up was done)
zmax  <- 5  # Lowest Depth Modelled
base  <- 1.05  # resolution reduction ("how strong is the grid corsened with depth")
# 1: equal spacing, >1 growing
dzmin <- 0.05 # dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]
Ti <- -15
omega=2*pi/(24*365*3600)
Seasonal_Amp=4
MAGST = -.5
#Difference between MAGST and Ti 14.5

#Experiment 2
dt <- 3600*24*5    #time step [s]
ns <- 1 #180  # how many time steps to compute before saving?
st <- 36500 # 365 # length of run [days]
zmax  <- 5  # depth of lowermost node center ("how large is the domain") [m]
base  <- 1.05  # resolution reduction ("how strong is the grid corsened with depth")
# 1: equal spacing, >1 growing
dzmin <- 0.05 # dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]
Ti <- -15
omega=2*pi/(24*365*3600)#missing 3600 in this formulation because of time step
Seasonal_Amp=4
MAGST = -5
#Difference between MAGST and Ti is 10 degrees

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
#tt=seq(1,nt, by=1)
##bc.upper <- rep(-.5, nt)+.08*t/(3600*24*365) 
bc.upper <- rep(MAGST, nt)+Seasonal_Amp*cos(omega*t)
##bc.upper <- rep(-.5, nt)+cos
  #seq(from = -1, to = 3, length.out = nt) #linear warming
bc.lower <- rep(0.2, nt)

#Material Properties
mat <- data.frame( zj = -soild$z, # depth of node [m]
                   dz = soild$dz, # with of node [m]
                   Tj = rep( Ti, nz), # initial temperature [C]
                   soi = rep(.5, nz), # volumetric proportion of soil matrix [0-1]
                   wat = rep(.5, nz), # volumetric proportion of total water [0-1]
                   liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
                   ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]
                   k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
                   c.soi = rep(1.8e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
)

#== Target arrays to hold results =========  this is where you run from after spinup
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
    out.t   <- c(out.t,t[i])
    print(paste(i,nt))
  }
}
## End of spin-up
##Beginning definition of plotting

plot(out.t/(3600*24*365),out.Tj[-731,19], main="Temperature", xlab="Years", 
     ylab="Temperature in Celsius",type="l")
plotdepth <- function(temp,mat,plotindex,time) {
  depth = as.character(round(mat$zj[plotindex],3))
  depth.nam = paste(depth,"m")
  plot(time/(3600*24*365),temp[-length(temp[,1]), plotindex],main=paste("ground temperature at ",depth.nam),
       xlab="years",
       ylab = "Temperature (C)", type='l')
}
#################
## Plot Depth, plots temperature vs time for specified depth index   
### usage: plotdepth(temp,mat,index,time)
##  where temp is a matrix of temperatures (usually out.Tj)
##  mat is a soil matrix ( could be mat or another)
## index is the 'number' of the depth 
##  time is a list of the times from the numerical model run

plotdepth(out.Tj,mat,3,out.t)
plot(out.t/(3600*24*365),out.Tj[-length(out.Tj[,5]),37], main="Temperature", xlab="Temperature in Celcius", 
     ylab="Depth z in m",type="l",ylim=c(-20 ,5))

for (i in c(1:length(out.Tj[,19]))){
  plot(out.Tj[i,],mat$zj, main="Temperature", xlab="Temperature in Celcius", 
       ylab="Depth z in m",type="l",xlim=c(-15 ,3))
  Sys.sleep(.1)
}

#This is where the changing of settings for warming begins
#== Settings for 2nd step of Part B:2 with Warming
dt <- 3600*24*5    #time step [s]
ns <- 1 #180  # how many time steps to compute before saving?
st <- 36500 # 365 # length of run [days]
# dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]


# initial temperature [C]

#== Preparation ================
st    <- st * 3600 * 24 # convert to [s]
t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
nt    <- length(t)
soild <- SoilDiscretize(dzmin, zmax, base)
nz    <- length(soild$z) #number of soil discretizations

backup.mat = mat
#make boundary condition
#bc.upper <- rep(-1.0, nt)
#tt=seq(1,nt, by=1)
##bc.upper <- rep(-.5, nt)+.08*t/(3600*24*365)
bc.upper <- rep(MAGST, nt)+Seasonal_Amp*cos(omega*t)+.08*t/(3600*24*365)
##bc.upper <- rep(-.5, nt)+cos
#seq(from = -1, to = 3, length.out = nt) #linear warming
bc.lower <- rep(0.2, nt)
#== Target arrays to hold results =========  this is where you run from after spinup
out.Tj.w  <- rbind(NULL,t(mat$Tj))
out.liq.w <- rbind(NULL,t(mat$liq))
out.t.w   <- t[0]
out.dz.w  <- mat$dz

#== Loop over time and solve =======
for (i in 1:nt) {
  mat <- HcGroundImplicit(dt, mat, bc.lower[i], bc.upper[i], unfrozen.type = unfrozen.type,
                          unfrozen.par = unfrozen.par, bcutype = bcutype, bcltype = bcltype,
                          layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
  #output only in specified interval
  check <- i/ns - floor(i/ns)
  if (check == 0) {
    out.Tj.w  <- rbind( out.Tj.w,t(mat$Tj))
    out.liq.w <- rbind(out.liq.w,t(mat$liq))
    out.t.w   <- c(out.t.w,t[i])
    print(paste(i,nt))
  }
}
##end of warming part
## Finding out when the permafrost is gone.

out.Tj.w[out.Tj.w[,1]>0,1] 
which(out.Tj.w[,1]>0)  # which gives you the index of a certain value

## what's the earliest index that temperature is above zero?
min(which(out.Tj.w[,1]>0))
## what about a little deeper?
min(which(out.Tj.w[,10]>0))
## what about a little deeper?
min(which(out.Tj.w[,20]>0))  ##  445 time steps
plotdepth(out.Tj.w,mat,20,out.t.w)  ## does it make sense?c
455*5/365  # how many years is that?

## now check all the depths!
min(which(out.Tj.w[,37]>0))  ##  445 time steps

# =============================================================================
#'
#' @title Numerical solution of ground heat conduction 
#'
#' @description Implicitly solves ground heat conduction based on finite-
#'              difference Crank-Nicolson scheme. Before the solution, ground
#'              and snow properties are updated  
#'
#' @details Positive are to the right while left shifts are expressed as a
#'          negative number. All shifts are circular. Elements shifted off one
#'          end wrap around and are shifted onto the other end. This function
#'          mimicks the behaviour of SHIFT in IDL. 
#'
#' @param dt Time step [s].
#' @param mat Data frame with relevant ground properties.  
#'
#' @param bc.lower Lower boundary condition as either [W m-2] (Neumann) or [C or K] 
#'            (Dirichlet). Array with dimension (N).  
#' @param bc.upper Upper boundary condition as either [W m-2] (Neumann) or [C or K] 
#'            (Dirichlet). Array with dimension (N). 
#' @param bcutype Type of upper boundary condition used. Can be "NEUMANN" 
#'                (standard, fixed heat flux) or "DIRICHLET" (fixed temperature).  
#' @param bcltype Type of upper boundary condition used. Can be "NEUMANN"
#'                (standard, fixed heat flux) or "DIRICHLET" (fixed temperature). 
#' 
#' @param unfrozen.type Character string indicating the type of function to be
#'                      used; the standard is "INTERVAL". The full options are: 
#"
#'                      "DALLAMICO": Dall’Amico, M., Endrizzi, S., Gruber, S., 
#'                      & Rigon, R. (2011). A robust and energy-conserving 
#'                      model of freezing variably-saturated soil. The 
#'                      Cryosphere, 5(2), 469–484. doi:10.5194/tc-5-469-2011
#'                      Typical values of van Genuchten parameters are found 
#'                      in Table 2 of Gubler, S., Endrizzi, S., Gruber, S., 
#'                      & Purves, R. S. (2013). Sensitivities and uncertainties 
#'                      of modeled ground temperatures in mountain environments. 
#'                      Geoscientific Model Development, 6(4), 1319–1336. 
#'                      doi:10.5194/gmd-6-1319-2013
#'
#'                      "MOTTAGHY": Mottaghy, D., & Rath, V. (2006). Latent heat 
#'                      effects in subsurface heat transport modelling and their
#'                      impact on palaeotemperature reconstructions. Geophysical 
#'                      Journal International, 164(1), 236-245.
#'
#'                      "INTERVAL": Phase change takes place in an interval of
#'                      specified with (unfrozen.par) below 0C and at a 
#'                      constant rate.  	
#'
#' @param unfrozen.par Parameter set for the chosen unfrozen water function.
#'
#'                     "DALLAMICO": unfrozen.par[1]: van Genuchten alpha [mm-1], 
#'                     unfrozen.par[2]: van Genuchten n [-], unfrozen.par[3]: 
#'                     residual water content [m3/m3]. The saturated water 
#'                     content is given by the input (mat$wat) and can thus 
#'                     vary with depth. 
#'
#'                     "MOTTAGHY": unfrozen.par[1]: width of freezing 
#'                     interval [K], unfrozen.par[2]: omega.  
#'
#'                     "INTERVAL" unfrozen.par[1]: width of freezing 
#'                     interval [K].  
#' @param type.sno Character string indicating the type of parameterization to be
#'             used for snow; the standard is CALONNE. Available are:
#'
#'             CALONNE: Calonne, N., Flin, F., Morin, S., Lesaffre, B., du 
#'             Roscoat, S. R., & Geindreau, C. (2011). Numerical and experimental
#'             investigations of the effective thermal conductivity of snow. 
#'             Geophysical Research Letters, 38(23), doi:10.1029/2011GL049234  
#'             "YEN": Yen, Y.-C. (1981). Review of thermal properties of 
#'              snow, ice and sea ice (34 pages). Hanover, NH, USA.
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             STURM: Sturm, M., J. Holmgren, M. König, and K. Morris (1997),
#'             The thermal conductivity of seasonal snow, Journal of 
#'             Glaciology, 43(143), 26–41.
#'
#'             JORDAN: Jordan, R. E., Andreas, E. L., & Makshtas, A. P. 
#'             (1999). Heat budget of snow-covered sea ice at North Pole 4. 
#'             Journal of Geophysical Research, 104(C4), 7785. 
#'             doi:10.1029/1999JC900011
#'
#'             WILLIAMS: Figure 4.11 in "The Frozen Earth: fundamentals of 
#'             geocryology" by P. Williams and M. Smith, 1989.
#' 
#' @param type.gnd Character string indicating the type of parameterization to be
#'             used for ground; the standard is "COSENZA". Available are:
#'
#'             COSENZA: Cosenza, P., Guerin, R., & Tabbagh, A. (2003). 
#'             Relationship between thermal conductivity and water content of 
#'             soils using numerical modelling. European Journal of Soil 
#'             Science, 54(3), 581–588. doi:10.1046/j.1365-2389.2003.00539.x
#'
#'             GEOMETRIC: Geometric mean, intermediate mixed conductivity model, 
#'             approximation of randomly oriented consituent elements. 
#'
#'             ARITHMETIC: Arithmetic mean, high mixed conductivity model, 
#'             approximation of consituent elements layered parallel to 
#'             temperature gradient.
#'
#'             HARMONIC: Harmonic mean, low mixed conductivity model, 
#'             approximation of consituent elements layered normal to 
#'             temperature gradient.
#'
#' @param layers.sno Index to identify snow layers within 'mat'.
#
#' @return Returns input array with Tj solved for next time step and temperature-
#'         dependent properties updated.
#' 
#' @export
#' @examples
#' #== Settings ==================
#' dt <- 3600    #time step [s]
#' ns <-  240    # how many time steps to compute before saving? 
#' st <- 4 * 365 # length of run [days]
#' dzmin <- 0.02 # minimal z-spacing ("how fine is the grid") [m]
#' zmax  <- 25   # depth of lowermost node center ("how large is the domain") [m]
#' base  <- 1.5  # resolution reduction ("how strong is the grid corsened with depth")
#'               # 1: equal spacing, >1 growing 
#' Ti <- -1      # initial temperature [C]
#'
#' type.gnd      <- "COSENZA"  # parameterization for soil thermal conductivity 
#' unfrozen.type <- "DALLAMICO" # paramterization for unfrozen water content
#' unfrozen.par  <- c(0.001, 1.4, 0.05)  # parameter(s) for unfrozen water content parameterization
#' bcutype       <- "DIRICHLET"# upper boundary condition type
#' bcltype       <- "NEUMANN"  # lower boundary condition type
#' layers.sno    <- c(0)       # index of nodes that are snow
#' type.sno      <- "CALONNE"  # parameterization for snow thermal conductivity
#'
#' #== Preparation ================
#' st    <- st * 3600 * 24 # convert to [s]
#' t     <- seq(from = 0, to = st, by = dt) #vector of time [s]
#' nt    <- length(t)
#' soild <- SoilDiscretize(dzmin, zmax, base)
#' nz    <- length(soild$z) #number of soil discretizations
#'
#' #make boundary condition
#' bc.upper <- seq(from = -1, to = 3, length.out = nt) #linear warming
#' bc.lower <- rep(0.05, nt)
#'
#' #make material properties and discretization
#' mat <- data.frame( zj = -soild$z, # depth of node [m]
#'                   dz = soild$dz, # with of node [m]
#'                   Tj = rep( Ti, nz), # initial temperature [C]
#'                  soi = rep(0.5, nz), # volumetric proportion of soil matrix [0-1]
#'                  wat = rep(0.0, nz), # volumetric proportion of total water [0-1]
#'                  liq = rep(  0, nz), # calculate later, volumetric proportion of liquid water [0-1]
#'                  ice = rep(  0, nz), # calculate later, volumetric proportion of ice [0-1]                
#'                k.soi = rep(2.5, nz), # thermal conductivity of soil matrix [W m-1 K-1]
#'                c.soi = rep(2e6, nz)  # volumetric heat capacity of soil matrix [J m-3 K-1]
#')
#'
#' #== Target arrays to hold results =========
#' out.Tj  <- rbind(NULL,t(mat$Tj))
#' out.liq <- rbind(NULL,t(mat$liq))
#' out.t   <- t[0]
#' out.dz  <- mat$dz
#'
#' #== Loop over time and solve =======
#' for (i in 1:nt) {
#'   mat <- HcGroundImplicit(dt, mat, bc.lower[i], bc.upper[i], unfrozen.type = unfrozen.type, 
#'                           unfrozen.par = unfrozen.par, bcutype = bcutype, bcltype = bcltype,
#'                           layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
                          
#' #output only in specified interval
#' check <- i/ns - floor(i/ns)
#'   if (check == 0) {
#'     out.Tj  <- rbind( out.Tj,t(mat$Tj))
#'     out.liq <- rbind(out.liq,t(mat$liq))
#'     out.t   <- t[i]
#'     print(paste(i,nt))   
#'   }		
#' }
#'
#' #== Plotting ================     
#' for (n in 1:147) {
#'	 plot(out.Tj[n,],mat$zj, type="l",xlim=c(-1.2 , 3))     
#'	 abline(v=0)
#' } 
#' 
#' @author Stephan Gruber <stephan.gruber@@carleton.ca>
#'
# =============================================================================

HcGroundImplicit <- function(dt, mat, bc.lower, bc.upper, 
                             unfrozen.type = "INTERVAL", unfrozen.par = c(1),
                             bcutype = "NEUMANN", bcltype = "NEUMANN",
                             layers.sno = c(0), type.sno = "CALONNE", type.gnd = "COSENZA") {
                             	
  #TODO: introduce heat transfer coefficient for boundary condition
  #make zero volumetric heat source term
  gj <- mat$Tj * 0
  
  #compute ground properties
  mat <- Unfrozen(mat, unfrozen.type, unfrozen.par)
  mat <- ThermCond(mat, layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
  mat <- HeatCapacity(mat)
  #mat <- SnowProp(mat, snow.type, snow.par)
  
  #TIME LOOP
  for (j in 1:length(bc.upper)) {

 	  
    #compute new solution (transpose to get from array to data frame)  	
  	mat$Tj <- Cranknicolson(mat$Tj, mat$kj, mat$cj, gj, mat$zj, dt, 
  	                          bc.lower[j], bc.upper[j], 
  	                          bcutype = bcutype, bcltype = bcltype)[,] 
  	#compute ground properties
    mat <- Unfrozen(mat, unfrozen.type, unfrozen.par)
    mat <- ThermCond(mat, layers.sno = layers.sno, type.sno = type.sno, type.gnd = type.gnd)
    mat <- HeatCapacity(mat)
    #mat <- SnowProp(mat, snow.type, snow.par)                                                                                               
  } 	  	  	  	
  return(mat)
}

# Analytical-exercise.R
#Original code in collaboration with Robert Martin
#Input Parameters - Make each of these inputs to a function later
Geotherm <- 0.2  
Vol_HC <- 1.8e6 #Volumetric Heat Capacity
Kt <- 2.5 #Thermal Conductivity
Period <- 365*24*3600
omega <- 2*pi/Period #How quickly the amplitude decreases w depth
Amp <- 4 
MAGST <- -0.5
kdiff <- Kt/Vol_HC #How easily the heat can move through the ground thermal diffusivity

n_time_steps=700
n_depth_steps=40

dt_an=5*3600*24
t_an=n_time_steps*dt_an

dz=.5
z=n_depth_steps*dz

#create the initial delT matrix  ## NOTE:  needs a mat$zj from numerical
delT=rep(0, length(mat$zj)*n_time_steps)
dim(delT)=c(length(mat$zj),n_time_steps)

for (z_int in c(1:length(mat$zj))){   ## need to run numerical first to get mat$zj
  #z=(z_int-1)*dz
  z = -mat$zj[z_int]
  Initial_Offset=MAGST + (Geotherm/Kt)*z
  for (t_int in 1:n_time_steps){
    
    t_an=(t_int-1)*dt_an
    
    delT[z_int,t_int] <- Initial_Offset +Amp * exp(-z*sqrt(omega/(2*kdiff)))*cos(omega*t_an-z*sqrt(omega/(2*kdiff)))
  }
}



plot((1:t_int)/73,delT[1,],type="l",xlab='years',ylim=c(-10,5))

lines((1:t_int)/73,delT[19,],type="l",xlab='years',col='red')


for (depthnumber in 2:n_depth_steps){
  lines((1:t_int)/73,delT[depthnumber,],type="l") #365/5= 73, 5 is the # of days in each timestep
  
  
}

#trumpet curve
#find the min, mean and average
#can pass a function the matrix with rows


