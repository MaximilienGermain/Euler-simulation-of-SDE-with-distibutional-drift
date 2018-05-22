% Computation of an approximation of a realisation of the solution of the 
% SDE dX = b(X)dt + dWt
function X = eulerMethod(X0,startNT,NT,N,T,H,B,xgrid,test,K,graphHaar)

% Setting the seed to 1
rng(1,'twister');

% Variables initialisation
startn = 1+ceil(T)*2^startNT; % Number of points on the first time grid
n = 1+ceil(T)*2^NT; % Number of points on the most refined time grid
dt = 2^(-NT); % Precision of the final grid
alpha = NT/startNT;

% Computation of the coefficients on the wavelet basis
Mu = computeMu(B,N,test,K);

% Display of the fBm and its derivative approximation
if (graphHaar == 1)
    figure 
    dispgrid=linspace(-K,K,3000);  
    approx = b(Mu,N,dispgrid);    
    [AX,H1,H2] = plotyy(dispgrid,approx,dispgrid,continuous(B,xgrid,dispgrid));
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
[t,W] = createfBm(0.5,ceil(T),NT-1,startNT-1,1+ceil(T)*2^startNT,0,2);
n=imax(T,t)+1;
t=t(1:n);
X = zeros(1,n);
X(1)=X0;
Xeuler = zeros(1,n);
if (test ~= 0)
    Xeuler(1) = X0;
end

for i=1:n-1
    X(i+1) = X(i) + b(Mu,N,X(i))*dt + W(i+1)-W(i);
    %X(i+1) = X(i) + b(Mu,N,X(i))*dt;
    if (test ~= 0)
        Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt + W(i+1)-W(i);
        %Xeuler(i+1) = Xeuler(i) + Xeuler(i)*dt;
    end
end

% Display
if (test == 0)
    figure
    modifiedEuler = plot(t,X,'b');
    grid on
    grid minor
    xlim([0 T])
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
    title(chn,'Interpreter','latex')
    
%     figure
%     modifiedEuler = plot(t,X-W(1:n)','b') ; %%%%% Test if there isn't saturation
%     grid on
%     grid minor
%         xlim([0 T])
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$X_t-W_t$','Interpreter','latex')
%     chn = ['Approximation of the drift part of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
% %   chn = ['Approximation of the drift part of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),' ; K = ',num2str(N),')'];
%     title(chn,'Interpreter','latex')

else
    figure
    modifiedEuler = plot(t,X,'r') ;
    grid on
    grid minor
    hold on
    xlim([0 T])
    euler = plot(t,Xeuler,'b');
    xlabel('$t$','Interpreter','latex')
    ylabel('$X_t$','Interpreter','latex')
    chn = ['Approximation of a sample path of the SDE solution (NT = ',num2str(NT),' ; N = ',num2str(N),')'];
    title(chn,'Interpreter','latex')
    legend([modifiedEuler,euler],'Approximation with Haar wavelets','Usual Euler scheme')
end

end