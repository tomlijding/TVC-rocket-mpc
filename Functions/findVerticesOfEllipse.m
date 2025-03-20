function vertices = findVerticesOfEllipse(P,c)
%findVerticesOfEllipse Finds the vertices of an Ellipse defined by x'*P*x =
%c
n = size(P,1);
bounds = zeros(n,1);
Pinv = P^-1;
for i = 1:n
    bounds(i) = sqrt(c*Pinv(i,i));
end
lb = -bounds;
ub = bounds;

% Number of dimensions (length of lb or ub)
numDims = length(lb);

% Generate all combinations of the bounds (binary choice for each dimension: lower or upper)
combinations = dec2bin(0:2^numDims-1) - '0';  % Generate binary combinations

% Initialize the vertices matrix
vertices = zeros(size(combinations, 1), numDims);

% Fill in the matrix by using the combinations of lower and upper bounds
for i = 1:numDims
    vertices(:, i) = lb(i) * (1 - combinations(:, i)) + ub(i) * combinations(:, i);
end

end

