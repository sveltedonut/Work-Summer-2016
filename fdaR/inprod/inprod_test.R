library(fda)
source("/Users/harry/Documents/MATLAB/fdaR-Harry/inprod/inprodC.R")
source("/Users/harry/Documents/MATLAB/fdaR-Harry/inprod/inprodT.R")

testbasis1 <- create.fourier.basis(c(0,1), 5)
testbasis2 <- create.fourier.basis(c(0,1), 5)

sq1 <- inprodT(testbasis1, testbasis2)
print(sq1)

sq2 <- inprodC(testbasis1, testbasis2)
print(sq1)

print(sq1 - sq2)