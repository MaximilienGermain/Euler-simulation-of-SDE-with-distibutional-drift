% Computes the coefficients of the B derivative approximation if test=0
% else it computes the coefficients of identity
function Mu = computeMu(B,N,test,K) %% Problem with N=1 ?????

Mu = zeros(N+2,K*2^(N+1)); 

if (test ==0)
    for m=1:2*K
        Mu(1,m) = -1/sqrt(2)*(B(getk(0,m-K-1,N,K))-B(getk(0,m-K,N,K)));
    end
    
    for j=2:N+2
        for m=1:K*2^(j-1)
            Mu(j,m) = -2^(j-2)*(B(getk(j-2,m-K*2^(j-2)-1,N,K))-2*B(getk(j-1,2*(m-K*2^(j-2)-1)+1,N,K))+B(getk(j-2,m-K*2^(j-2)-1+1,N,K)));
        end
    end
else
    for m = 1:2*K
        Mu(1,m) = (2*(m-K-1)+1)/(2*sqrt(2));
    end    
    
    for j=2:N+2
        Mu(j,:) = - 1/(2^(j));
    end
    
end

end