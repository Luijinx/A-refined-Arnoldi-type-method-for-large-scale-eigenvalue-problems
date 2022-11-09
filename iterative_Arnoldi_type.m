function [e, V, res]=iterative_Arnoldi_type(A, m)
n=2000;
%A=sprand(n, n, 6/n);
%A=A+A';
q=rand(n, 1); q=q/norm(q);
i=0;
k=4; tau=6; 

tol=1e-6;
iter=70;
residuo=-1;

while(i<=iter && residuo==-1)
    i=i+1;
    
    [e, V, res]=Arnoldi_type(A, m, k, tau, q);

    for j=1:k
        if(res(j)>=tol)
            residuo=-1;
        else 
            residuo=0;
        end
    end
    
    %formo un nuovo vettore q
    if(residuo==-1)
        q=zeros(n, 1);
        for j=1:k
        q=q+res(j)*imag(V(:, j))+real(V(:, j));
        end
    end
end
