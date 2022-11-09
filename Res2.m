function [a]=Res2(A, G, Q)
  n=size(A, 1);
  [V, E]=eig(G'*G);
  [~, p]=min(diag(E));
  a=V(:, p);
  a=Q*a;
  a=a/norm(a);
end