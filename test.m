% To fix omega, we must fix seed and startNT. To fix b we must fix startN, and Kmax
X0 = 1;
startNT = 10;
startN = 3; 
NT = 10; % Time grid precision 2^(-N)
N = 8; % Space grid precision 2^(-N)
T = 1.2;
Kmax = 5;
Nx = 1+Kmax*2^(startN+2); % 2 times more precise than the grid for b^N
H = 0.85;
graphHaar = 0;
control = 0;
testId = 0;
seed = 2;
PlotActive = 1;

[xgrid,B,~] = createfBm(H,Kmax,N,startN,Nx,-Kmax,1000);
Mu = computeMu(B,N,testId,Kmax);
rng('shuffle');

graphHaar = 1;
% control = 0;
[X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive);
imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.eps'])

% Test with b=id
% testId = 1;
% graphHaar = 1;
% Kmax = 5;
% N = 4;
% B = linspace(-Kmax,Kmax,1+Kmax*2^(N+2));
% Mu = computeMu(B.^2/2,N,testId,Kmax);
% [X,~,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B.^2/2,Mu,B,testId,Kmax,graphHaar,control,seed,PlotActive);
% imshow(['Haar H = ',num2str(H),' ; N = ',num2str(N),'.png'])
% imshow(['Identity test n = ',num2str(NT),' ; N = ',num2str(N),'.png'])

% Convergence in n
% startNT = 2;
% PlotActive = 1;
% %control = 1;
% graphHaar = 1;
% nmin = 2;
% nmax = 12;
% convergencen(nmin,nmax,X0,startNT,N,T,H,B,Mu,xgrid,testId,Kmax,graphHaar,control,seed,PlotActive)

% Convergence in N
% control = 1;
% graphHaar = 1;
% PlotActive = 1;
% minN = startN;
% maxN = startN+5;
% convergence_N(minN,maxN,X0,startNT,startN,NT,T,H,testId,Kmax,seed,Nx,graphHaar,control,PlotActive)

%  Monte-Carlo N
% MC = 1000;
% NT = 12; % Error depends a lot on NT 
% startNT = NT;
% Nmax = 9;
% Kmax = 6;
% minN = Nmax-5;
% startN = minN;  
% maxN = Nmax-1; %%   
% PlotActive = 0;
% [expectations,var] = monteCarlo_N(X0,T,H,Nmax,Kmax,graphHaar,control,testId,minN,maxN,PlotActive,MC,NT,startNT,startN)

%  Monte-Carlo n
% M = 400;
% X0 = 0;
% NTmax = 12;
% minN = 5;
% maxN = NTmax-2;
% PlotActive = 0;
% [expectations,var] = monteCarlon(X0,xgrid,B,N,T,Mu,H,Kmax,graphHaar,control,testId,minN,maxN,PlotActive,M,NTmax)