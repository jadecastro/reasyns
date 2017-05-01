function [ex] = sampleEllipsoidBoundary(c,Q,d)

if length(c) ~= length(Q)
    error('Dimensions of c and Q must be the same!')
end

n = length(c);

[V,E] = eig(Q);
e = sqrt(diag(E));

% Normalize the ellipsoid to have a smallest axis of 1.
ai = e./min(e);

ex = [];
while (isempty(ex))
    % Generate points with uniform density on a sphere surface
    u = rand(d,n-1);
    theta(:,1:n-2) = 2*pi*u(:,1:n-2);
    theta(:,n-1) = acos(2*u(:,n-1)-1);
    theta(1:d,n) = 0;
    
    if n == 2
        X(:,1) = sin(theta(:,1)).*cos(theta(:,2));
        X(:,2) = cos(theta(:,1));
    elseif n == 3
        X(:,1) = sin(theta(:,1)).*sin(theta(:,2)).*cos(theta(:,3));
        X(:,2) = sin(theta(:,1)).*cos(theta(:,2));
        X(:,3) = cos(theta(:,1));
    elseif n == 4
        X(:,1) = sin(theta(:,1)).*sin(theta(:,2)).*sin(theta(:,3)).*cos(theta(:,4));
        X(:,2) = sin(theta(:,1)).*sin(theta(:,2)).*cos(theta(:,3));
        X(:,3) = sin(theta(:,1)).*cos(theta(:,2));
        X(:,4) = cos(theta(:,1));
    elseif n == 5
        X(:,1) = sin(theta(:,1)).*sin(theta(:,2)).*sin(theta(:,3)).*sin(theta(:,4)).*cos(theta(:,5));
        X(:,1) = sin(theta(:,1)).*sin(theta(:,2)).*sin(theta(:,3)).*cos(theta(:,4));
        X(:,2) = sin(theta(:,1)).*sin(theta(:,2)).*cos(theta(:,3));
        X(:,3) = sin(theta(:,1)).*cos(theta(:,2));
        X(:,4) = cos(theta(:,1));
    else
        error('n = 1 or n > 5 not implemented.')
    end
    % Select points with probability inversely related to how far they are from
    % the surface of the original sphere. The acceptance rate appears to be
    % about 50% in the worst case of an elongated cigar shape.
    t = rand(d,1);
    term = X.^2./repmat(ai',d,1);
    idx = find((sum(term,2) >= t.^2));
    acceptance_rate = 100*size(idx)/size(t);
    
    % stretch the selected points onto the ellipse
    ex = X(idx,:).*repmat(e',length(idx),1);
    ex = (V*ex')' + repmat(c',length(idx),1);
end
