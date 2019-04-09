% Created by Maximiliano Rivera, based on Numerical Linear Algebra by
% L.N.Trefetehn.

function [ Q, R ] = QRbyGSM(A)
% This function does QR decomposition with Gram-Schmidt modify algorithm

[m, n] = size(A);
Q = zeros(m ,n);
R = zeros(n, n);
for i=1:n
    vi = A(:, i);
    rii = norm(vi);
    qi = vi/rii;
    Q(:, i) = qi;
%     for j=i+1:n
%         vj = A(:, j);
%         rij = qi'*vj;
%         A(:, j) = vj - rij*qi;
%     end
    Vj = A(:, i+1:end);
    rij = qi'*Vj;
    A(:, i+1:end) = Vj - qi*rij;
    
    % We can construct R from the rii and rij previously calculated
    R(i:n, i) = [rii; rij'];
end
R = R';
return

end