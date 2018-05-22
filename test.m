% To fix omega, we must fix startNT, start, and Kmax
X0 = 1;
startNT = 1;
NT = 8; % Time grid precision 2^(-N)
N = 5; % Space grid precision 2^(-N)
start = 1; 
T = 1.2;
%K = ceil(exp(2)+T*3);
%K = start+5;
Kmax = start+8;
Nx = 1+Kmax*2^(start+2); % 2 times more precise than the grid for b^N
H = 0.55;
[xgrid,B,M] = createfBm(H,Kmax,N,start,Nx,-Kmax,1000);

X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,1);

% Test with b=id
% B = linspace(-Kmax,Kmax,1+Kmax*2^(N+2));
% X = eulerMethod(X0,startNT,NT,N,T,H,B.^2/2,B,1,Kmax,1);

% Convergence in n
% for NT=6:10
%     X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,0);
% end

% Convergence in N
% for N=3:7
%     [xgrid,B,M] = createfBm(H,K,7,start,Nx,-K,1000);
%     NT = 8;
%     X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,0,Kmax,0);
% end