function [ Q, R ] = QRbyGSC(A)
% This function is for decompose a matrix A in matrices Q and R, where Q
% has orthonormal columns and R is upper triangular, such that A = QR. This is possible
% using the Gran-Schmidt Clasical algorithm

[m, n] = size(A);
Q = zeros(m ,n);
R = zeros(n, n);

for j=1:n
    v = A(:, j);
    r = conj(Q)'*v;
    v = v - Q*r;
    rjj = norm(v);
    v = v/rjj;
    Q(:, j) = v; 
    R(1:j, j) = [r(1:j-1); rjj];
    
    
end
return

end