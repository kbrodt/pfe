function [Ka, Ma, N]=EFMaxwell_P1(S,T,NAr,na,pos,mu_1,mu_2, epsilon_1,epsilon_2)
nt=size(T,1); ns=size(S,1);
Ka=sparse(na,na);
Ma=sparse(na,na);
N = sparse(2*ns,na);
for t=1:nt
    [Kat, Mat, Nt]=EFelmMaxwell_P1(S(T(t,:),:), sign(NAr(t,:)), pos(t),mu_1,mu_2, epsilon_1,epsilon_2);
    I=T(t,:);K=abs(NAr(t,:));
    Ka(K,K)=Ka(K,K)+Kat;
    Ma(K,K)=Ma(K,K)+Mat;
    N([I I+ns],K)=N([I I+ns],K)+Nt;
end

function [Ka, Ma, N]=EFelmMaxwell_P1(S,sA,pos,mu_1,mu_2, epsilon_1,epsilon_2)
nbq=7;os=sqrt(15);s3=1./3.;
pp1=(6.-os)/21.;pp2=(6.+os)/21.;pp3=(9.+2.*os)/21.;pp4=(9.-2.*os)/21.;
pts_quadT=[s3 s3; pp1 pp1; pp1 pp3; pp3 pp1;pp2 pp2; pp2 pp4; pp4 pp2];
pp1=(155.-os)/2400.;pp2=(155.+os)/2400.;
pds_quadT=[9./80.;pp1;pp1;pp1;pp2;pp2;pp2];
S21=S(2,:)-S(1,:);S31=S(3,:)-S(1,:);
delta=S21(1)*S31(2)-S21(2)*S31(1);aire=abs(delta)/2;
Jflmt=[S31(2) -S21(2);-S31(1) S21(1)]/delta;
Ma=zeros(3,3);
N=zeros(6,3);
rotr=2*sA/delta; Ka=1./mu(pos,mu_1,mu_2)*aire*rotr'*rotr;
for k=1:nbq,
    x=pts_quadT(k,1);y=pts_quadT(k,2);
    pk=pds_quadT(k)*abs(delta);
    P=S(1:3,:)'*[1-x-y;x;y];
    w=[1-x-y x y];
    x=P(1);y=P(2);
    r=[sA;sA].*[S(3,2)-y S(1,2)-y S(2,2)-y;x-S(3,1) x-S(1,1) x-S(2,1)]/delta;
    Ma=Ma+epsilon(pos,epsilon_1,epsilon_2)*pk*r'*r;
    N=N+[pk*w'*r(1,:); pk*w'*r(2,:)];
end,