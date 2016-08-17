function [AB]=aretes_bord(A,NA,B)
na=size(A,1);nt=size(NA,1);AB=zeros(na,1);
nb=size(B,1);
for a=1:na
    for b=1:nb
        if (sum( A(a,:)==B(b,:) )==2)
            AB(a)=1;
        end
    end,
end