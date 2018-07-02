% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function [X,Y,frame,haar,control,usual] = eulerMethod(X0,startNT,NT,N,T,H,B,Mu,xgrid,testId,K,graphHaar,control,seed,plot)

frame = [];
haar = [];
control = [];
usual = [];

% Precision of the time grid
dt = 2^(-NT); 

% Brownian motion simulation on the time grid
kmax=floor(2^startNT*T);
[t,W] = createfBm(0.5,(kmax+1)/2^(startNT),NT-1,startNT-1,kmax+2,0,seed);
n = length(t);

% Initialisation of X (and Xeuler)
X = zeros(1,n);
X(1) = X0;
Y = zeros(1,n);
Y(1) = X0;
Xeuler = zeros(1,n);
if (testId ~= 0)
    Xeuler(1) = X0;
end
lastX = X0;
lastY = X0;

% Euler scheme
for i=1:n-1
    if abs(lastX) > K || abs(lastY) > K
        t = t(1:i);
        X = X(1:i);
        Y = Y(1:i);
        break
    end
    inc = W(i+1) - W(i);
    lastX = lastX + b(Mu,K,lastX)*dt + inc;
    lastY = lastY + b(Mu,K,lastY)*dt - inc;
    X(i+1) = lastX;
    Y(i+1) = lastY;
    if (testId ~= 0)
        Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + W(i+1)-W(i);
    end
end

% Displays of the approximated solution
if (plot ~= 0)
    if (testId == 0)
        frame = plotSolution(t,X,NT,N,T);
        %plotSolution(t,Y,NT,N,T);
    else
        usual = plotUsualSolution(t,X,NT,N,T,Xeuler);
    end
end

if (control == 1)
    control = plotControl(t,X,NT,N,T,W);
end

% Display of the fBm and its derivative approximation
if (graphHaar == 1)
    haar = plotHaarApproximation(K,Mu,N,B,xgrid,H);
end

end