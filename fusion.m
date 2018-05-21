function c = fusion(a,b)

la=length(a);
c=zeros(la*2,1);
for i=1:la
    c(2*i-1) = a(i);
    c(2*i) = b(i);
end

end