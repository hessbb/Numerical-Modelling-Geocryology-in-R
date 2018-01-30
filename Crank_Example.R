Tj <- array(c(  1,    1,    1,    1), dim = c(1, 4))
kj <- array(c(  2,    2,    2,    2), dim = c(1, 4))
cj <- array(c(2e6,  2e6,  2e6,  2e6), dim = c(1, 4))
gj <- array(c(  0,    0,    0,    0), dim = c(1, 4))
zj <- array(c(  0,   -1,   -2,   -3), dim = c(1, 4))
dt <- 30000
bcl <-  3
bcu <- -4

#set up plot
plot(Tj, zj, type="l", xlim=c(-5,5))

#loop over heat conduction and add new lines
for (n in 1:100) {
  Tj <- Cranknicolson(Tj, kj, cj, gj, zj, dt, bcl, bcu,
                      bcutype = "DIRICHLET", bcltype = "NEUMANN")
  lines(Tj, zj)
}