% Haar basis functions h_{j,m} (x)
function output = h(j,m,x)

% Test if j is a integer greater than -1
if (~(j == floor(j))) || (j < -1)
    error('j must be an integer greater than -1 when you calcuate h(j,m,x)')
end

y = x*2^j - m;
output = zeros(1,length(x));

for k=1:length(x)
    
    % Test if y is inside the haar function support
    if (y(k) >= 1) || (y(k) < 0)
        output(k) = 0;
    else
        if (j == -1) || (y(k) < 1/2)
            output(k) = 1;
        else   
            if (y(k) >= 1/2) 
                output(k) = -1; 
            end
        end   
    end
    
end

end