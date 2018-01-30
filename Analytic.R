#Input Parameters - Make each of these inputs to a function later
Analytic <- function(Vol_HC=1.8 *1000000,Kt= 2.5,Period= 365*24*3600,Amp= 4){
  

  
    



omega <- 2*pi/Period
kdiff <- Kt/Vol_HC

dt=10*3600*24
t=71*dt

dz=.25
z=100*dz

#create the initial delT matrix
delT=rep(0, 100*71)
dim(delT)=c(100,71)

for (z_int in 1:100){
  z=z_int*dz
  
  for (t_int in 1:71){
    
    t=t_int*dt
    
    delT[z_int,t_int] <- Amp * exp(-z*sqrt(omega/(2*kdiff)))*cos(omega*t-z*sqrt(omega/(2*kdiff)))
  
  }
  
}



plot(1:t_int,delT[1,],type="l")

for (pl in 2:100){
  lines(1:t_int,delT[pl,],type="l")
  
  #trumpet curve
  #find the min, mean and average
  #can pass a function the matrix with rows  
}
}




