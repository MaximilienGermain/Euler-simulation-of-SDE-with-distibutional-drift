function y=continuous(B,xgrid,x)

n=length(x);
N=length(xgrid);
y=zeros(n,1);

for i=1:n
    if (x(i)>=xgrid(1))||(x(i)<=xgrid(N))
        for j=1:N-1
            if (x(i)>= xgrid(j))&&(x(i)<xgrid(j+1))
                y(i) = B(j)+(x(i)- xgrid(j))/(xgrid(j+1)-xgrid(j))*(B(j+1)-B(j));
            end
        end
    end
end

end