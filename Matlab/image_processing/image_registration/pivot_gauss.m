%%
% CS 6640 : Image Processing Project 3
%
% Author : Arthur COSTE
% Date : October 2012
%
% Content : Gaussian elimination
%% 
function x = pivot_gauss(A,b)

%rank of the matrix
n=size(A,1);

for p=1:n
    vec=[(1:p-1) n (p:n-1)];
    test=1;
    %test if matrix is invertible
    while A(p,p)==0
        if test==n
            error('Cannot invert matrix')
        end
        A=A(vec,:);
        b=b(vec);
        test=test+1;
    end
    %perform division by the pivot
    b(p)=b(p)/A(p,p);
    A(p,:)=A(p,:)/A(p,p);
    % perform substraction of rows
    for q=p+1:n
        b(q)=b(q)-A(q,p)*b(p);
        A(q,:)=A(q,:)-A(q,p)*A(p,:);
    end
end
%compute values for unknowns
x=zeros(n,1);
x(n)=b(n);
%from bottom to top
for p=n-1:-1:1
    x(p)=b(p);
    for q=p+1:n
    x(p)=x(p)-A(p,q)*x(q);
    end
end