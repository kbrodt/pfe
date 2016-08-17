function [AA,FF]=elimine(AA,FF,Refneu)
N=length(AA);
for i=1:N
    if Refneu(i)~=0
        Ad=AA(i,i);
        AA(i,:)=0;
        AA(:,i)=0;
        AA(i,i)=Ad;
        FF(i)=0;
    end
end