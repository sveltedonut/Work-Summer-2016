addpath('./')
addpath('../../ProfIRT')
addpath('../../fdaM')

load thetaH
load WfdH
load UH
load INDH
load TQH
load SlambdaH
load CritH
load ItermaxH

dataStr.U   = UH;
dataStr.ind = INDH;
dataStr.tQ  = TQH;

mex thetafunC.c bsplineC.c

tic;
theta1      = thetafun(thetaH, WfdH, dataStr, SlambdaH);
toc

tic;
theta2      = thetaC(thetaH, WfdH, dataStr, SlambdaH);
toc

max(max(abs(theta1 - theta2)))