% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function [X,frame] = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,testId,K,graphHaar,control)

% Setting the seed to 1 for reuse
rng(1,'twister');

% Precision of the time grid
dt = 2^(-NT); 

% Computation of the coefficients on the Haar wavelet basis
Mu = computeMu(B,N,testId,K);

% Brownian motion simulation on the time grid
kmax=floor(2^startNT*T);
[t,W] = createfBm(0.5,(kmax+1)/2^(startNT),NT-1,startNT-1,kmax+2,0,2);
n = length(t);

% Initialisation of X (and Xeuler)
X = zeros(1,n);
X(1) = X0;
Xeuler = zeros(1,n);
if (testId ~= 0)
    Xeuler(1) = X0;
end

% Euler scheme
for i=1:n-1
    X(i+1) = X(i) + b(Mu,N,X(i))*dt + W(i+1)-W(i);
    if (testId ~= 0)
        Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + W(i+1)-W(i);
    end
end

% Displays of the approximated solution
if (testId == 0)
    frame = plotSolution(t,X,NT,N,T)
else
    plotUsualSolution(t,X,NT,N,T,Xeuler)
end

if (control == 1)
    plotControl(t,X,NT,N,T,W)
end

% Display of the fBm and its derivative approximation
if (graphHaar == 1)
    plotHaarApproximation(K,Mu,N,B,xgrid,H)
end

end