addpath('/Users/harry/Documents/MATLAB/fdaM')
addpath('/Users/harry/Documents/MATLAB/fdaM-Harry/bspline')

testbasis   = create_bspline_basis([0,2], 100, 4);
evalarg     = 0:0.01:2;

%tic;
basismat1   = getbasismatrix(evalarg, testbasis, 0)
%toc

%tic;
basismatC   = getbasismatrixC(evalarg, testbasis, 0)
%toc

subplot(2, 1, 1);
plot(basismat1);
subplot(2, 1, 2);
plot(basismatC);
max(max(abs(basismat1 - basismatC)))