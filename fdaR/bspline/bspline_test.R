library(fda)
source("/Users/harry/Documents/MATLAB/fdaR-Harry/bspline/getbasismatrixC.R")

testbasis <- create.bspline.basis(c(0,1), 5)
evalarg <- seq(0, 1, 0.1)

basismat1 <- getbasismatrix(evalarg, testbasis, 0)
print(basismat1)

basismat2 <- getbasismatrixC(evalarg, testbasis, 0)
print(basismat1)

print(basismat1 - basismat2)