function [a]=Res(S, lambda, Vm)
  n=size(S, 1);
  W=(S-lambda*speye(n))*Vm;
  [V, E]=eig(W'*W);
  [~, p]=min(diag(E));
  a=V(:, p);
  a=Vm*a;
  a=a/norm(a);
end
