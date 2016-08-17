function [B1el, B2el, B3el, B4el] = mat_elem_surface(S1, S2, ref)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul les matrices de masse surfaciques pour les frontieres 
% Gamma1, Gamma2, Gamma3 et Gamma4
%
% SYNOPSIS [B1el, B2el, B3el, B4el] = mat_elem_surface(S1, S2, ref)
%          
% INPUT * S1, S2 : les 2 coordonnees des 2 sommets de l'arete 
%                      (vecteurs reels 1x2)
%         ref    : Reference de l'arete.
%
% OUTPUT - B1el matrice de masse surfacique elementaire pour Gamma1 (matrice 2x2)
%		 - B2el matrice de masse surfacique elementaire pour Gamma2 (matrice 2x2)
%        - B3el matrice de masse surfacique elementaire pour Gamma3 (matrice 2x2)
%        - B4el matrice de masse surfacique elementaire pour Gamma4 (matrice 2x2)
%
% NOTE  le calcul est exacte (pas de condensation de masse).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);

% On calcul la longueur de l'arete.
Long = sqrt((x2-x1)^2+(y2-y1)^2);

% Declaration et ramplisage des matrice elementaires.
%---------------------------------------------------
B1el = zeros(2,2);
B2el = zeros(2,2);
B3el = zeros(2,2);
B4el = zeros(2,2);

if ref==1      % On est sur Gamma1
	% A COMPLETER
    B1el = [1/3 1/6; 1/6 1/3]*Long;
elseif ref==2      % On est sur Gamma2
	% A COMPLETER
    B2el = [1/3 1/6; 1/6 1/3]*Long;
elseif ref==3  % On est sur Gamma3
	% A COMPLETER
    B3el = [1/3 1/6; 1/6 1/3]*Long;
elseif ref==4  % On est sur Gamma4
	% A COMPLETER
    B4el = [1/3 1/6; 1/6 1/3]*Long;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
