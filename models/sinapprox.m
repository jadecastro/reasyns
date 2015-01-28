function y = sinapprox(x,n)

y = 0;
% x = mod(x+pi,2*pi) - pi;

for i = 1:n
    pm = -(-1)^i;
    y = y + pm*x^(2*i-1)/factorial(2*i-1);
end
