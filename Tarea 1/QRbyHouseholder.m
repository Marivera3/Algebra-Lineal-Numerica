% Created by Maximiliano Rivera, based on Numerical Linear Algebra by
% L.N.Trefetehn.

function [ Q, R ] = QRbyHouseholder(A)

[m, n] = size(A);
Q = eye(m);
for k = 1:n
   
    x = A(k:m, k);
    e1 = zeros(length(x), 1);
    e1(1) = 1;
    v = sign(x(1))*norm(x)*e1 + x;
    v = v/norm(v);
    F = eye(m-k+1, m-k+1) - 2*v*(v');
    
    A(k:m, k:n) = F*A(k:m, k:n);
    
    Qk = eye(m, m);
    Qk(end-(m-k):end, end-(m-k):end) = F; 
    Q = Qk*Q; % forming the matrix Q' = Qk...Q2Q1
end
Q = Q';
R = A;
return
end