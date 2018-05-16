% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function X = eulerMethod(X0,NT,N,T,Mu)

% Setting the seed to 1
rng(1,'twister');

% Variables initialisation
dt = T/NT;
t = T/NT*(0:NT);
X = zeros(1,NT+1);
Xeuler = zeros(1,NT+1);
X(1) = X0;
Xeuler(1) = X0;

for i=1:NT
    increment = sqrt(dt)*randn();
    X(i+1) = X(i) + b(Mu,X(i))*dt + increment;
    Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + increment;
end

% Display
figure
modifiedEuler = plot(t,X,'--or') ;
hold on
euler = plot(t,Xeuler,'--ob');
xlabel('t')
ylabel('X_t')
chn = ['Approximation of one realisation of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
title(chn)
legend([modifiedEuler,euler],'Approximation with Haar basis','Usual Euler scheme')

end