function [Q, R, U]=Arnoldi_2(A, m, q, tau)
n=size(A, 1);
q=q/norm(q);
U=zeros(n, m);
Q=zeros(n, m);
R=zeros(m, m);
G=zeros(m, m-1);

U(:, 1)=q;
R(1, 1)=norm((A-tau*eye(n))*q); Q(:, 1)=(A-tau*eye(n))*q/R(1, 1);

for j=1:m-1
    z=(A-tau*eye(n))*Q(:, j);
    for i=1:j
        G(i, j)=dot(z, Q(:, i));
        z=z-G(i, j)*Q(:, i);
    end
    
    G(j+1, j)=norm(z);
    if(G(j+1, j)<=1e-10)
        disp('K_j è invariante')
        return
     end
    
     
    Q(:, j+1)=z/G(j+1, j);
end

R(:, 2:m)=G;
U(:, 2:m)=Q(:, 1:m-1);

        