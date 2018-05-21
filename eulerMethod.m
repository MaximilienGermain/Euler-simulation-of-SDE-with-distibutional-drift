% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,test,K)

% Setting the seed to 1
rng(1,'twister');

% Variables initialisation
startn = 1+startNT; % Number of points on the first time grid
n = 1+NT; % Number of points on the most refined time grid
dt = T/NT; % Precision of the final grid
t = linspace(0,T,n);
X = zeros(1,n);
Xeuler = zeros(1,n);
X(1) = X0;
if (test ~= 0)
    Xeuler(1) = X0;
end

% Computation of the coefficients on the wavelet basis
Mu = computeMu(B,N,test,K);

% Display of the fBm and its derivative approximation
% if (test == 0)
%     figure
%     approx = b(Mu,K,xgrid);
%     [AX,H1,H2] = plotyy(xgrid,approx,xgrid,B);
%     ylabel(AX(1),'Approximation of $\frac{\partial}{\partial x}B^H_x$ by Haar wavelets','Interpreter','latex')
%     ylabel(AX(2),'$B^H_x$','Interpreter','latex')
%     set(AX(1),'ycolor','r')
%     set(AX(2),'ycolor','b')
%     set(H1,'Color','r')
%     set(H2,'Color','b')
%     set(AX(2),'Xgrid','on','Ygrid','on','YMinorGrid','on')
%     xlabel('$x$','Interpreter','latex')
%     xlim(AX(2), [min(xgrid) max(xgrid)])
%     xlim(AX(1), [min(xgrid) max(xgrid)])
%     ylim(AX(1), [min(approx)-0.2*abs(min(approx)) max(approx)+0.2*abs(min(approx))])
%     ylim(AX(2), [min(B)-0.2*abs(min(B)) max(B)+0.2*abs(min(B))])
%     chn = ['Sample path of $B^H_x$ and approximation of $\frac{\partial}{\partial x}B^H_x$ (H = ',num2str(H),' ; N = ',num2str(N),')'];
%     title(chn,'Interpreter','latex')
% end

% Euler scheme
rnd = randn(1,startn-1);
%for p=1:NT-startNT
for p=1:NT/startNT-1
    newrandom = randn(1,2^(p-1)*(startn-1));
    temprandom = zeros(1,2^(p)*(startn-1));
    for y=1:2^(p-1)*(startn-1)
        temprandom(2*y-1) = newrandom(y);
        temprandom(2*y) = rnd(y);
    end
    rnd = temprandom;
end
rnd = sqrt(dt)*rnd;
%n = 1+(T)*2^NT;

W = [0,cumsum(rnd)];

for i=1:n-1
    X(i+1) = X(i) + b(Mu,K,X(i))*dt + rnd(i);
    %X(i+1) = X(i) + b(Mu,K,X(i))*dt;
    if (test ~= 0)
        %Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + rnd(i);
        Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt;
    end
end

% Display
% if (test == 0)
%     figure
%     modifiedEuler = plot(t,X,'b');
%     grid on
%     grid minor
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$X_t$','Interpreter','latex')
%     chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),' ; K = ',num2str(K),')'];
%     title(chn,'Interpreter','latex')
%     %legend([modifiedEuler],'Approximation with Haar wavelets')
%     
%     figure
%     modifiedEuler = plot(t,X-W,'b') ; %%%%% Test if there isn't saturation
% else
%     figure
%     modifiedEuler = plot(t,X,'r') ;
%     grid on
%     grid minor
%     hold on
%     euler = plot(t,Xeuler,'b');
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$X_t$','Interpreter','latex')
%     chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),' ; K = ',num2str(K),')'];
%     title(chn,'Interpreter','latex')
%     legend([modifiedEuler,euler],'Approximation with Haar wavelets','Usual Euler scheme')
% end

end