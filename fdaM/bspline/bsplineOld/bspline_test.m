addpath('/Users/harry/Documents/MATLAB/fdaM')
addpath('/Users/harry/Documents/MATLAB/fdaM-Harry/bspline')

testbasis   = create_bspline_basis([0,2], 4, 2);
testbasisL  = create_bspline_basis([0,2], 4, 1);
evalarg     = 0:0.1:2;

%tic;
basismat1   = getbasismatrix(evalarg, testbasis, 1)


basismatL   = getbasismatrixC(evalarg, testbasisL, 0)
%toc
%tic;
basismat0   = getbasismatrix(evalarg, testbasis, 0)
%toc

%tic;
basismatC   = getbasismatrixC(evalarg, testbasis, 1)
%toc

subplot(2, 2, 1);
plot(basismat1);
subplot(2, 2, 3);
plot(basismat0);
subplot(2, 2, 4);
plot(basismatL);
subplot(2, 2, 2);
plot(basismatC);
max(max(abs(basismat1 - basismatC)))