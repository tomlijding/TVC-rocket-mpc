function outputMatrix = jitter(inputMatrix)
%JITTER Function makes matrices positive semi-definite after numerical
%errors.
N = size(inputMatrix,2);
jitterMatrix = 1e-6*eye(N);
outputMatrix = inputMatrix + jitterMatrix;
end

