function y = tanapprox(x,n)

% Note: abs(x) < pi/2

y = 0;

for i = 1:n
    B = bernoullinumber(2*i);
    y = y + B*(-4)^i*(1-4^i)*x^(2*i-1)/factorial(2*i);
end
