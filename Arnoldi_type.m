function [e, V, res]=Arnoldi_type(A, m, k, tau, q)
%m non può essere troppo alto altrimenti è nu burdell

%inizializzazione elementi
n = size(A,1);
res=zeros(1, k);
%q=rand(n, 1);
q=q/norm(q);

%creazione delle matrici Q, R, U, B
[Q, R, U]=Arnoldi_2(A, m, q, tau);

B=R/((Q')*U);


%risolvo il sistema con la matrice B
[V, D] = eig(B);

e=diag(D);


[val, ind]=sort(abs(e));
e=e(ind);
V=V(:, ind);


%scelta dei k autovalori tra m (prendo i più piccoli perchè sono quelli
%più vicini al punto tau)

e=e(1:k);
V=V(:, 1:k);
e=e+tau;
V=Q*V;


%calcolo del residuo
for j=1:k
     res(j) = norm(A*V(:,j) - e(j)*V(:,j))/norm(A, 1);
end




