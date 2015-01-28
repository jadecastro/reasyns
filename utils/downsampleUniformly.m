function [Y,indices] = downsampleUniformly(X,N)
% Downsample a monotonic vector X approximately uniformly based on a
% desired factor N (i.e. keep every N sample). Return both the downsampled
% vector Y and the vector of sampled indices.

X0 = X(1);
X = X - X0;
NN = length(X);
delX = (N/NN)*(X(end) - X(1));
indices(1) = 1;
for i = 1:NN/N
    if any(X > i*delX)
        indices(i+1) = find(X >= i*delX,1,'first');
    end
end
indices(end+1) = NN;
indices = unique(indices);
Y = X(indices) + X0;
