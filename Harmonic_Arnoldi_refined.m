function [mv, E, V, res]=Harmonic_Arnoldi_refined(A, m, k, tau, v1)
%function [cpu, mv, E, V, res]=Harmonic_Arnoldi_refined(A, m, k, tau, v1)
% Approssima i k autovalori della matrice A piu' vicini a tau mediante
% Gli autovalori della proiezione di A nello spazio di Krylov di dimensione m-1
% span(v1, A v1, A^2 v1,... A^(m-2)v1) 

%In questa versione cambia solo la maniera di calcolare gli autovettori
   
   n = size(A,1);
   E = zeros(k, 1);
   V = zeros(n, k);
%Computiamo le matrici H tilde, H_m e V_m con il metodo di Arnoldi,
%calcolando anche il numero mv 
   [~, H, Vm, Hm, ~, ~, mv] = Arnoldi(A, m, v1);

%Procedo utilizzando Rayleigh-Ritz per computare psi e theta
   B1=(Hm-tau*eye(m))';
   I = eye(m+1); I = I(1:m+1, 1:m);
   B2 = (H-tau*I);
   B2 = B2'*B2;
   e = eig(B1, B2);
% scelgo i k autovalori di modulo più grande e calcolo col metodo di
% minimizzazione gli autovettori corrispondenti

   [val,ind] = sort(abs(e),'descend');
   e = e(ind);
   E = e(1:k).^(-1)+tau;
   for j=1:k
      V(:, j)=Res(A, E(j), Vm);
% check res
      w = A*V(:,j); 
      res(j) = norm(w - E(j)*V(:,j));
   end
