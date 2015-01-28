function y = cosapprox(x,n)

y = 1;
% x = mod(x+pi,2*pi) - pi;

for i = 1:n
    pm = (-1)^i;
    y = y + pm*x^(2*i)/factorial(2*i);
end
