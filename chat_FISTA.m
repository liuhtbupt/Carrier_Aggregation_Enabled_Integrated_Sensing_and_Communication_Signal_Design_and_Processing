% FISTA algorithm for compressed sparse data reconstruction
% Input: A: measurement matrix
%        b: measurement data
%        lambda: regularization parameter
% Output: x: reconstructed sparse signal

function x = chat_FISTA(A,b,lambda)

% Set the maximum number of iterations and the tolerance
max_iter = 1000;
tol = 1e-6;

% Convert A to sparse matrix if it is sparse
if nnz(A)/numel(A) < 0.1 % check the sparsity ratio of A
    A = sparse(A); % convert A to sparse matrix
end

% Initialize
[m,n] = size(A); % get the size of A
x = zeros(n,1); % initialize x as a zero vector
y = x; % initialize y as x
t = 1; % initialize t as 1
L = norm(A)^2; % calculate the Lipschitz constant of A

% Main loop
for k = 1:max_iter
    
    % Gradient descent step
    z = y - (A'*(A*y-b))/L;
    
    % Soft thresholding step
    x_new = sign(z).*max(abs(z)-lambda/L,0);
    
    % Update t and y
    t_new = (1+sqrt(1+4*t^2))/2;
    y = x_new + (t-1)/t_new*(x_new-x);
    
    % Check the stopping criterion
    if norm(x_new-x)/norm(x) < tol
        break;
    end
    
    % Update x and t
    x = x_new;
    t = t_new;
    
end

end
