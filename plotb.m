x = 1/100*(-700:700);
N = 4;

% Construction of the coefficients matrix
Mu = zeros(N+2,K*2^(N+1));
for m = 1:2*K+1
    Mu(1,m) = (2*(m-K-1)+1)/(2*sqrt(2));
end    
    
for j=2:N+2
    Mu(j,:) = - 1/(2^(j));
end

value = b(Mu,K,x);
figure
plot(x,value)
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$Id^N(x)$','Interpreter','latex')
chn = ['Approximation of the identity function by Haar wavelets (N = ',num2str(N),')'];
title(chn,'Interpreter','latex')