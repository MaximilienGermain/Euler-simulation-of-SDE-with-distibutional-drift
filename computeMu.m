% Computes the coefficients of the B derivative approximation if test=0
% else it computes the coefficients of identity
function [Mu,Muerr] = computeMu(B,N,test,Kmax,H) %% Problem with N=1 ?????

Mu = zeros(N+2,Kmax*2^(N+1)); 
Muerr = zeros(N+2,Kmax*2^(N+1)); 

for m=1:2*Kmax
    Mu(1,m) = -(B(getk(0,m-Kmax-1,N,Kmax))-B(getk(0,m-Kmax,N,Kmax)));
    Muerr(1,m) = -2^(1-H+1/2)*(B(getk(0,m-Kmax-1,N,Kmax))-B(getk(0,m-Kmax,N,Kmax)));
end
    
for j=2:N+2
    for m=1:Kmax*2^(j-1)
        Mu(j,m) = -(B(getk(j-2,m-Kmax*2^(j-2)-1,N,Kmax))-2*B(getk(j-1,2*(m-Kmax*2^(j-2)-1)+1,N,Kmax))+B(getk(j-2,m-Kmax*2^(j-2)-1+1,N,Kmax)));
        Muerr(j,m) = -2^((j-2)*(H-1/2))*(B(getk(j-2,m-Kmax*2^(j-2)-1,N,Kmax))-2*B(getk(j-1,2*(m-Kmax*2^(j-2)-1)+1,N,Kmax))+B(getk(j-2,m-Kmax*2^(j-2)-1+1,N,Kmax)));
        %Mu(j,m) = -2^(j-2)*(B(getk(j-2,m-Kmax*2^(j-2)-1,N,Kmax))-2*B(getk(j-1,2*(m-Kmax*2^(j-2)-1)+1,N,Kmax))+B(getk(j-2,m-Kmax*2^(j-2)-1+1,N,Kmax)));
    end
end
end