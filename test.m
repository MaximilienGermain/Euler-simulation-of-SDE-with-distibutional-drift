% To fix omega, we must fix startNT, start, and Kmax
X0 = 1;
startNT = 1;
NT = 6; % Time grid precision 2^(-N)
N = 5; % Space grid precision 2^(-N)
start = 1; 
T = 1.2;
Kmax = start+9;
Nx = 1+Kmax*2^(start+2); % 2 times more precise than the grid for b^N
H = 0.85;
graphHaar = 0;
control = 0;
[xgrid,B,M] = createfBm(H,Kmax,N,start,Nx,-Kmax,1000);

%
% graphHaar = 0;
% control = 1;
% X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,graphHaar,control);

% Test with b=id
% graphHaar = 1;
% B = linspace(-Kmax,Kmax,1+Kmax*2^(N+2));
% X = eulerMethod(X0,startNT,NT,N,T,H,B.^2/2,B,0,Kmax,graphHaar,control);

% Convergence in n
% for NT=6:10
%     startNT = 6;
%     X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,0,control);
% end

% Convergence in N
for N=1:10
    [xgrid,B,M] = createfBm(H,Kmax,7,start,Nx,-Kmax,1000);
    NT = 8;
    control = 1;
    X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,graphHaar,control);
end