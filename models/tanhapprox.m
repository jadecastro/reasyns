function y = tanhapprox(x,n)

y = 0;

for i = 1:n
    B = bernoullinumber(2*i);
    y = y + B*4^i*(4^i-1)*x.^(2*i-1)/factorial(2*i);
end
