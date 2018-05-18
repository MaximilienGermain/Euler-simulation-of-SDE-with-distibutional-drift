% Here newM is the the transpose of the cholesky root
function [newrandom,newgrid,newB,newM] = refinefbm(random,xgrid,H,N,M)

rng(100,'twister');
Nx = length(xgrid);
newgrid = linspace(min(xgrid),max(xgrid),2*Nx-1);
dx = 1/2^(N+2);
otherM = M';

% Simulation of Nx - 1 independent gaussian variables
otherrandom = randn(1,Nx-1)';

for u=1:Nx-1
   s = length(otherM);
   alpha = zeros(1,s);
   for k=1:s+1-u
       alpha(k) = dx^(2*H)*((2*u-1)^(2*H)+(2*k)^(2*H)-abs(2*u-1-2*k)^(2*H))/2;
   end
   for k=s+2-u:s % the real k goes from 1 to u-1 here
       alpha(k) = dx^(2*H)*((2*u-1)^(2*H)+(2*(k-s+u-1)-1)^(2*H)-abs(2*u-1-2*(k-s+u-1)+1)^(2*H))/2;
   end
   %V = 100*alpha/(100*otherM)';
   V = alpha/(otherM)';
   w = sqrt(dx^(2*H)*(2*u-1)^(2*H)-V*V');
   otherM = [otherM zeros(s,1); V w];
end
gamma = otherM*otherM';
allrandom = [random ; otherrandom];
newvalues = otherM*allrandom;
newB = [0];
%newB = zeros(1,2*length(newvalues)+1);

newgamma = zeros(2*Nx-2);
for i=1:Nx-1
    newgamma(2*i-1,:) = gamma(Nx-1+i,:);
    newgamma(2*i,:) = gamma(i,:);    
end
gamma = newgamma;
for j=1:Nx-1
    newgamma(:,2*j-1) = gamma(:,Nx-1+j);
    newgamma(:,2*j) = gamma(:,j);    
end
newM = chol(newgamma);

newrandom = zeros(2*Nx-2,1);
for t=1:Nx-1
    newrandom(2*t-1) = random(t);
    newrandom(2*t) = otherrandom(t);
end

for o=1:length(newvalues)/2
    newB(2*o) = newvalues(o+Nx-1);
    newB(2*o+1) = newvalues(o);
end
% Sort newM maybe in order to use refinefbm several times ?

end
