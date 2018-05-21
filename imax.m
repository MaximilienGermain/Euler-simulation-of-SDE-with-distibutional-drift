function im = imax(T,t)

n=length(t);
for i=1:n
   if (t(i)<=T) && (t(i+1)>T)
       im=i;
   end
end

end