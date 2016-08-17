function [A,NA]=aretes(S,T)
ns=size(S,1); nt=size(T,1); na=0;
G=sparse(ns,ns); NA=zeros(nt,3);
for t=1:nt,
    for a=1:3,
        ap=a+1; if(ap==4) ap=1;end,
        if(T(t,a)<T(t,ap)) I=T(t,a);J=T(t,ap);s=1;
        else I=T(t,ap); J=T(t,a);s=-1;
        end
        I=min(T(t,a),T(t,ap)); J=max(T(t,a),T(t,ap));
        if(G(I,J)==0)
            na=na+1; G(I,J)=na;
        end,
        NA(t,a)=s*G(I,J);
    end,
end,
A=zeros(na,2);[r,c,v]=find(G);
A(abs(v),1)=r;A(abs(v),2)=c;