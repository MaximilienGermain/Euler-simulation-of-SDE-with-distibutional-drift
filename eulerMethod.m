% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function X = eulerMethod(X0,NT,N,Nx,T,H,B,xgrid,test)

% Setting the seed to 1
rng(1,'twister');

% Variables initialisation
dt = T/NT;
t = T/NT*(0:NT);
X = zeros(1,NT+1);
Xeuler = zeros(1,NT+1);
X(1) = X0;
Xeuler(1) = X0;

Mu = computeMu(B,N,test);
figure
[AX,H1,H2] = plotyy(xgrid,b(Mu,xgrid),xgrid,B);
ylabel(AX(1),'Approximation of $\frac{\partial}{\partial x}B^H_x$','Interpreter','latex')
ylabel(AX(2),'$B^H_x$','Interpreter','latex')
set(AX(1),'ycolor','r')
set(AX(2),'ycolor','b')
set(H1,'Color','r')
set(H2,'Color','b')
set(AX(2),'Xgrid','on','Ygrid','on','YMinorGrid','on')
xlabel('$x$','Interpreter','latex')
xlim([min(xgrid) max(xgrid)])
chn = ['Sample path of $B^H_x$ and approximation of its derivative (Nx = ',num2str(Nx),' ; H = ',num2str(H),')'];
title(chn,'Interpreter','latex')

for i=1:NT
    increment = sqrt(dt)*randn();
    X(i+1) = X(i) + b(Mu,X(i))*dt + increment;
    Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + increment;
end

% Display
if (test == 0)
    figure
    modifiedEuler = plot(t,X,'r') ;
    grid on
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
    title(chn,'Interpreter','latex')
    legend([modifiedEuler],'Approximation with Haar wavelets')
else
    figure
    modifiedEuler = plot(t,X,'r') ;
    grid on
    hold on
    euler = plot(t,Xeuler,'b');
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
    title(chn,'Interpreter','latex')
    legend([modifiedEuler,euler],'Approximation with Haar wavelets','Usual Euler scheme')
end