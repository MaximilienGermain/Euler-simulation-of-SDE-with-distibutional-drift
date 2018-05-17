% Here newM is the the transpose of the cholesky root
function [newgrid,newB,newM] = refinefbm(random,xgrid,H,N,M)

rng(1000,'twister');
Nx = length(xgrid);
newgrid = linspace(min(xgrid),max(xgrid),2*Nx-1); %%% To check
dx = 1/2^(N+2);
newM = M';

% Simulation of Nx - 1 independent gaussian variables
newrandom = randn(1,Nx-1)';

for u=1:Nx-1
   s = length(newM);
   alpha = zeros(1,s);
   for k=1:s+1-u
       alpha(k) = dx^(2*H)*((2*u-1)^(2*H)+(2*k)^(2*H)-abs(2*u-1-2*k)^(2*H))/2;
   end
   for k=s+2-u:s
       alpha(k) = dx^(2*H)*((2*u-1)^(2*H)+(2*(k-s+u-1)-1)^(2*H)-abs(2*u-1-2*(k-s+u-1)+1)^(2*H))/2;
   end
   V = alpha/newM';
   w = sqrt(dx^(2*H)*(2*u-1)^(2*H)-V*V');
   newM = [newM zeros(s,1); V w];
end

newvalues = newM*[random ; newrandom];
newB = [0];
for o=1:length(newvalues)/2
    newB(2*o) = newvalues(o+Nx-1);
    newB(2*o+1) = newvalues(o);
end

end
