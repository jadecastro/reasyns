function [df, d2f, d3f] = dynamicsGradientsDubins(a1, a2, a3, a4, order)
% This is an auto-generated file.
%
% See <a href="matlab: help generateGradients">generateGradients</a>. 

% Check inputs:
typecheck(a1,'DubinsPlant');
if (nargin<4) order=1; end
if (order<1) error(' order must be >= 1'); end
sizecheck(a1,[1  1]);
sizecheck(a2,[1  1]);
sizecheck(a3,[3  1]);
sizecheck(a4,[1  1]);

% Symbol table:
a3_3=a3(3);
a4_1=a4(1);


% Compute Gradients:
df = sparse(3,5);
df(1,4) = -(2*sin(a3_3))/25;
df(2,4) = (2*cos(a3_3))/25;
df(3,5) = 1;

% d2f
if (order>=2)
  d2f = sparse(3,25);
  d2f(1,19) = -(2*cos(a3_3))/25;
  d2f(2,19) = -(2*sin(a3_3))/25;
else
  d2f=[];
end  % if (order>=2)

% d3f
if (order>=3)
  d3f = sparse(3,125);
  d3f(1,94) = (2*sin(a3_3))/25;
  d3f(2,94) = -(2*cos(a3_3))/25;
else
  d3f=[];
end  % if (order>=3)


 % NOTEST
