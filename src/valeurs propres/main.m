% =====================================================
%
% ---------------------------------
nom_maillage = 'geometrie1.msh';

mu_1 = 1;
mu_2 = 1;
epsilon_1 = 1;
epsilon_2 = 0.01; % sigma = epsiolon_1^(-1) = 10;

% lecture du maillage et affichage
% ---------------------------------
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,NbaretesB,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KKn = sparse(Nbpt,Nbpt); % matrice de rigidite
MMn = sparse(Nbpt,Nbpt); % matrice de rigidite
FFn = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  
  % calcul des matrices elementaires du triangle l 
  
   [Kel]=matK_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:));
   % LA ROUTINE matK_elem.m DOIT ETRE MODIFIEE

   [Mel]=matM_elem(Coorneu(Numtri(l,1),:),Coorneu(Numtri(l,2),:),...
		       Coorneu(Numtri(l,3),:));
    
    % On fait l'assemblage de la matrice globale
    % A COMPLETER  
    I=Numtri(l,1);J=Numtri(l,2);K=Numtri(l,3);
    MMn([I,J,K],[I,J,K])=MMn([I,J,K],[I,J,K])+Mel;
    KKn([I,J,K],[I,J,K])=KKn([I,J,K],[I,J,K])+Kel; 
end % for l

% Matrice EF
% -------------------------

% pseude-elimination pour matrice de rigidite neoudes
[KKn,tilde_FF]=elimine(KKn,FFn,Refneu);
% figure
% spy(KKn);

[MMn,tilde_FF]=elimine(MMn,FFn,Refneu);
%figure
%spy(AA^-1);

% tableaux des aretes
[Ar,NumtriAr]=aretes(Coorneu,Numtri);

% Matrice de masse et rigidite des elements d'aretes 
na=size(Ar,1);
[KKa,MMa,N] = EFMaxwell_P1(Coorneu,Numtri,NumtriAr,na,Reftri,mu_1,mu_2,epsilon_1,epsilon_2); 
% figure
% spy(MMa);

% renumeroter des aretes sur le bord, car ils ne sont pas bon apres "lecture_msh"
for i=1:NbaretesB
    if(Numaretes(i,1)>Numaretes(i,2))
        tmp = Numaretes(i,1);
        Numaretes(i,1)=Numaretes(i,2);
        Numaretes(i,2)=tmp;
    end
end

% label 1  pour les aretes sur le bord 
[RefArBord]=aretes_bord(Ar,NumtriAr,Numaretes);

% pseude-elimination pour matrice de rigidite d'aretes
IdArBord=find(RefArBord==1); 
IdArInt = find(~RefArBord);
KKaInt = KKa(IdArInt,IdArInt);
% DA=diag(KKa);
% KKa(IdArBord,:)=0; KKa(:,IdArBord)=0;
% for i=1:size(IdArBord,1),KKa(IdArBord(i),IdArBord(i))=DA(IdArBord(i));end

% pseude-elimination pour matrice de masse d'aretes
MMaInt = MMa(IdArInt,IdArInt);
% DA=diag(MMa);
% MMa(IdArBord,:)=0; MMa(:,IdArBord)=0;
% for i=1:size(IdArBord,1),MMa(IdArBord(i),IdArBord(i))=DA(IdArBord(i));end
% figure
% spy(MMa);

% Resolution du pb aux v.p. KKa*E = l*MMa*E
% nbvp=300;
% [V0,D0]=eigs(KKa,MMa,nbvp,'sm');
% [lambda0,Ix0]=sort(diag(D0));

% operateur Gradient transpose ou bien un operateur, qui return le bord des aretes (summ de neoudes) : Aretes -> Neoudes
Gt = GradT(Nbpt,Ar);
% figure
% spy(Gt);

% нужно сделать псевдо-элиминацию дл€ оператора √рад“, то есть нужно
% обнулить все столбцы Gt, которые соответсвуют ребрам на границе <-
% не работает !!! нужно добавить и те ребра, которые имеют одну вершину на
% границу, а другую внутри области...
% pseudo-elimination d'operateru GradT
% pas de relation entre aretes de BORD et INTERIEUR du domaine


KKna = sparse(Nbpt,Nbpt);
KKna = Gt*MMa*Gt';

% IdArTouchBord=find(Ar(:,1)<=NbaretesB & Ar(:,2)>NbaretesB);
% % Gt(:,IdArTouchBord)=0; % PAS BON!!!
% for i=1:size(IdArTouchBord,1)
%   IdNeuBord=find(Gt(:,IdArTouchBord(i))==-1);
%   Gt(IdNeuBord,IdArTouchBord(i))=0;
% end
IdSomInt = find(~Refneu);
% Gt(:,IdArBord)=0;
% Gt(find(Refneu),:)=0;


KKnaInt = KKna(IdSomInt,IdSomInt);
% IdArInt = find(Ar(:,1)>NbaretesB & Ar(:,2)>NbaretesB);
% GtInt = Gt(IdSomInt,:);
% GinvKKnGtIn = GtInt'*KKnaInt^-1*GtInt;
% GinvKKnGtInt = GinvKKnGtIn(IdArInt,IdArInt);

GtInt = Gt(IdSomInt,IdArInt);
%GtInt = Gt(IdSomInt,:);
GinvKKnGtInt = GtInt'*KKnaInt^-1*GtInt;
%GinvKKnGtInt = GinvKKnGtInt(IdArInt,IdArInt);

% NInt=N([IdSomInt IdSomInt+Nbpt],IdArInt);
NInt=N(IdSomInt,IdArInt);

% [KKna,tilde_FF]=elimine(KKna,FFn,Refneu);
% figure
% spy(KKna);
% GinvKKnGt = Gt'*KKna^-1*Gt;
% DA=diag(GinvKKnGt);
% GinvKKnGt(IdArBord,:)=0; GinvKKnGt(:,IdArBord)=0;
% for i=1:size(IdArBord,1),GinvKKnGt(IdArBord(i),IdArBord(i))=DA(IdArBord(i));end
% 
% DA=diag(MMa);
% MMa(IdArBord,:)=0; MMa(:,IdArBord)=0;
% for i=1:size(IdArBord,1),MMa(IdArBord(i),IdArBord(i))=DA(IdArBord(i));end

% PPa = MMa-MMa*GinvKKnGt*MMa;
PPaInt = MMaInt-MMaInt*GinvKKnGtInt*MMaInt;
PP = eye(size(PPaInt,1))-GinvKKnGtInt*MMaInt;
% PPaInt = PPa(IdArInt,IdArInt);
% DA=diag(PPa);
% PPa(IdArBord,:)=0; PPa(:,IdArBord)=0;
% for i=1:size(IdArBord,1),PPa(IdArBord(i),IdArBord(i))=DA(IdArBord(i));end
% IdArInt = find(~RefArBord);
% % IdArInt = find(Ar(:,1)>NbaretesB & Ar(:,2)>NbaretesB);
% MMaInt = MMa(IdArInt,IdArInt);
% KKaInt = KKa(IdArInt,IdArInt);
% PPaInt = PPa(IdArInt,IdArInt);

% Resolution du pb aux v.p. KKa*E = l*PPa*E
nbvp=250;%size(IdArInt,1);
% MMaInt=(MMaInt+MMaInt')/2;
% KKaInt=(KKaInt+KKaInt')/2;
% PPaInt=(PPaInt+PPaInt')/2;%+eye(size(PPaInt,1))*1e-10;
% PPaInt=PPaInt+eye(size(PPaInt,1))*1e-10;
% [V,D]=eigs(KKaInt,PPaInt,nbvp,'sm');
% [lambda,Ix]=sort(diag(D));
% PPa*E=l*MMa*E
GAUCHE = KKaInt+MMaInt;%NInt'*KKnaInt^-1*NInt;%KKaInt+eye(size(PPaInt,1))*1e-10;%
DROITE = MMaInt;%+NInt'*KKnaInt^-1*GtInt;%MMaInt*PP;%MMaInt*GtInt'*KKnaInt^-1*NInt;
% GAUCHE=(GAUCHE+GAUCHE')/2;
% DROITE=(DROITE+DROITE')/2;
[VRight,D]=eig(full(GAUCHE),full(DROITE*PP));
% [WLeft,D]=eig(full(GAUCHE'),full(DROITE'));
[lambda,Ix]=sort(diag(D));
lambda=lambda-1;

PEprRight= PP*VRight;
normERight=zeros(size(PPaInt,1),1);
normPERight=zeros(size(PPaInt,1),1);
ErrRelRight=zeros(size(PPaInt,1),1);
% PEprLeft= PP*WLeft;
% normELeft=zeros(size(PPaInt,1),1);
% normPELeft=zeros(size(PPaInt,1),1);
% ErrRelLeft=zeros(size(PPaInt,1),1);
% cosTheta=zeros(size(PPaInt,1),1);
for ii=1:size(PPaInt,1)
    normPERight(ii)=norm(PEprRight(:,ii));
    normERight(ii)=norm(VRight(:,ii));
    ErrRelRight(ii) = norm(PEprRight(:,ii)-VRight(:,ii))/normERight(ii);
%     cosTheta(ii)=dot(VRight(:,ii),WLeft(:,ii))/(norm(VRight(:,ii))*norm(WLeft(:,ii)));
end
normPERight=normPERight(Ix);
normERight=normERight(Ix);
ErrRelRight=ErrRelRight(Ix);
% cosTheta=cosTheta(Ix);
% Residu = KKn - KKna;
% figure
% spy(Residu);

% n'iprote quoi ----------
% IdArB=find(Ar(:,1)<=NbaretesB & Ar(:,2)<=NbaretesB);
% Gt(:,IdArB)=0;
% I = [];
% for i=1:NbaretesB;
%     ArBord = Ar(IdArBord(i),:);
%     Nl = [];
%     Nr = [];
%     for j=1:size(IdArTouchBord,1)
%         ArBordT = Ar(IdArTouchBord(j),:);
%         if(ArBord(1)==ArBordT(1))
%             Nl = [Nl ArBordT(2)];
%         end
%         if(ArBord(2)==ArBordT(1))
%             Nr = [Nr ArBordT(2)];
%         end
%     end
%     I = [I intersect(Nl,Nr)];
% end
% for i=1:NbaretesB
%     Gt(I(i),IdArB(i))=1;
% end

%IdNeu=Ar(IdAr,2);
%Gt(IdNeu,IdAr)=1;
% figure
% spy(Gt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

