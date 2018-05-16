x = 1/100*(-500:500);
N = 15;

% Construction of the coefficients matrix
Mu = zeros(N+2,2*N+1);
for m = 1:2*N+1
    Mu(1,m) = (2*(m-N-1)+1)/(2*sqrt(2));
end    
    
for j=2:N+2
    Mu(j,:) = - 1/(2^(j));
end

value = b(Mu,x);
figure
plot(x,value)
grid on
xlabel('$x$','Interpreter','latex')
ylabel('$Id^N(x)$','Interpreter','latex')
chn = ['Approximation of the identity function by Haar wavelets (N = ',num2str(N),')'];
title(chn,'Interpreter','latex')