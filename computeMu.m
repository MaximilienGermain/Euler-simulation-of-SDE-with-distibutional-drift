function Mu = computeMu(B,N,test)

Mu = zeros(N+2,2*N+1); 

if (test ==0)
    for m=1:2*N+1
        Mu(1,m) = -1/sqrt(2)*(B(m)-B(m+1));
    end
    
    for j=2:N+2
        for m=1:2*N+1
            Mu(j,m) = -2^(j)*(B(getk(j-2,m-N,N))-2*B(getk(j-1,2*(m-N-1)+1,N))+B(getk(j-2,m-N-1+1,N)));
        end
    end
else
    for m=1:2*N+1
        Mu(1,m) = (2*(m-N-1)+1)/(2*sqrt(2));
    end    
    
    for j=2:N+2
        Mu(j,:) = - 1/(2^(j));
    end
    
end

end