function Gt=GradT(ns,A)
Gt = sparse(ns,size(A,1));
for j=1:size(A,1)
        Gt([A(j,1),A(j,2)],j) = [-1; 1];
end