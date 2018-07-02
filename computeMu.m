% Computes the coefficients of the B derivative approximation if test=0
% else it computes the coefficients of identity
function Mu = computeMu(B,N,test,Kmax) %% Problem with N=1 ?????

Mu = zeros(N+2,Kmax*2^(N+1)); 

%if (test ==0)
    for m=1:2*Kmax
        Mu(1,m) = -1/sqrt(2)*(B(getk(0,m-Kmax-1,N,Kmax))-B(getk(0,m-Kmax,N,Kmax)));
    end
    
    for j=2:N+2
        for m=1:Kmax*2^(j-1)
            Mu(j,m) = -2^(j-2)*(B(getk(j-2,m-Kmax*2^(j-2)-1,N,Kmax))-2*B(getk(j-1,2*(m-Kmax*2^(j-2)-1)+1,N,Kmax))+B(getk(j-2,m-Kmax*2^(j-2)-1+1,N,Kmax)));
        end
    end
% else
%     for m = 1:2*N
%         Mu(1,m) = (2*(m-N-1)+1)/(2*sqrt(2));
%     end    
%     
%     for j=2:N+2
%         Mu(j,:) = - 1/(2^(j));
%     end    
% end

end