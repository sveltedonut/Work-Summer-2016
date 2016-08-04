addpath('/Users/harry/Documents/MATLAB/fdaM')

testbasis1 = create_fourier_basis([0,1], 200);
testbasis2 = create_bspline_basis([0,1], 200);
%harmaccelLfd = harmaccel(testbasis);
%harmaccelLfd2 = harmaccel(testbasis2);

% plot(testbasis);

%testfdPar = fdPar(testbasis, harmaccelLfd, 1e-7);
%testLfd = getLfd(testfdPar);
%testfdPar2 = fdPar(testbasis2, harmaccelLfd2, 1e-7);
%testLfd2 = getLfd(testfdPar2);

tic;
sq1 = inprod(testbasis1, testbasis2);
toc

tic;
sq2 = inprodC(testbasis1, testbasis2);
toc

sq1 - sq2;
%plot(sq2);