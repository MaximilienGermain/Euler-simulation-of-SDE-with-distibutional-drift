% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function X = eulerMethod(X0,NT,N,T,H,B,xgrid,test,K)

% Setting the seed to 1
rng(1,'twister');

% Variables initialisation
dt = T/NT;
t = T/NT*(0:NT);
X = zeros(1,NT+1);
Xeuler = zeros(1,NT+1);
X(1) = X0;
if (test ~= 0)
    Xeuler(1) = X0;
end

% Computation of the coefficients on the wavelet basis
Mu = computeMu(B,N,test,K);

% Display of the fBm and its derivative approximation
if (test == 0)
    figure
    approx = b(Mu,K,xgrid);
    [AX,H1,H2] = plotyy(xgrid,approx,xgrid,B);
    ylabel(AX(1),'Approximation of $\frac{\partial}{\partial x}B^H_x$ by Haar wavelets','Interpreter','latex')
    ylabel(AX(2),'$B^H_x$','Interpreter','latex')
    set(AX(1),'ycolor','r')
    set(AX(2),'ycolor','b')
    set(H1,'Color','r')
    set(H2,'Color','b')
    set(AX(2),'Xgrid','on','Ygrid','on','YMinorGrid','on')
    xlabel('$x$','Interpreter','latex')
    xlim(AX(2), [min(xgrid) max(xgrid)])
    xlim(AX(1), [min(xgrid) max(xgrid)])
    ylim(AX(1), [min(approx)-0.2*abs(min(approx)) max(approx)+0.2*abs(min(approx))])
    ylim(AX(2), [min(B)-0.2*abs(min(B)) max(B)+0.2*abs(min(B))])
    chn = ['Sample path of $B^H_x$ and approximation of $\frac{\partial}{\partial x}B^H_x$ (H = ',num2str(H),' ; N = ',num2str(N),')'];
    title(chn,'Interpreter','latex')
end

% Euler scheme
for i=1:NT
    increment = sqrt(dt)*randn();
    X(i+1) = X(i) + b(Mu,K,X(i))*dt + increment;
    if (test ~= 0)
        Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + increment;
    end
end

% Display
if (test == 0)
    figure
    modifiedEuler = plot(t,X,'b') ;
    grid on
    grid minor
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),' ; K = ',num2str(K),')'];
    title(chn,'Interpreter','latex')
    %legend([modifiedEuler],'Approximation with Haar wavelets')
else
    figure
    modifiedEuler = plot(t,X,'r') ;
    grid on
    grid minor
    hold on
    euler = plot(t,Xeuler,'b');
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),' ; K = ',num2str(K),')'];
    title(chn,'Interpreter','latex')
    legend([modifiedEuler,euler],'Approximation with Haar wavelets','Usual Euler scheme')
end

end