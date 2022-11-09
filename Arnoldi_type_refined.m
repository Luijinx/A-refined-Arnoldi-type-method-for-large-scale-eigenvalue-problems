function [e, V, res]=Arnoldi_type_refined(A, m, k, tau, q)

%inizializzazione elementi
n = size(A,1);
V=zeros(n, k);
res=zeros(1, k);
%q=rand(n, 1); 
q=q/norm(q);
   
%creazione delle matrici Q, R, U, B
[Q, R, U]=Arnoldi_2(A, m, q, tau);
B=R/((Q')*U);

cond(
%risolvo il sistema con la matrice B
[V1, D] = eig(B);

e=diag(D);

[val, ind]=sort(abs(e));
e=e(ind);

%scelta dei k autovalori tra m
e=e(1:k);
e=e+tau;


%calcolo del residuo
for j=1:k
    G=U/R-Q/(e(j)-tau);
    V(:, j)=Res2(A, G, Q);
     res(j) = norm(A*V(:,j) - e(j)*V(:,j))/norm(A, 1);
end

