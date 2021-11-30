function [elementosOrd] = Elements_sorter(nodos,elementos)

% Funcion que ordena los nodos del elemento H8 según el criterio:
% ksi //= x
% eta //= y
% zeta//= z

            %%%%%%%%% Orientación del elemento %%%%%%%%%%%%
% 
%                                 zeta // (z) 
%                                  | 
%                             2 ---|------ 1
%                            /|    |      /          
%                           / |    x     /|
%                          /  |         / |
%                         3 ---------- 4 x-------- eta // (y)
%                         |   |        |  |   
%                         |   6 --x----|- 5
%                         |  /   /     | /
%                         | /   /      |/
%                         7 ---/------ 8
%                             /
%                           ksi // (x)

%% Inicio del ordenamiento:
elementosOrd = zeros(size(elementos));

for iele = 1:size(elementos,1);    
    
nodosEle = nodos(elementos(iele,:),:); 
%% Ordeno por Z:
[coord1,indiceZ] = sort(nodosEle(:,3));

IndexCaraArribaZ = indiceZ(5:8);
IndexCaraAbajoZ = indiceZ(1:4);

%% Para la cara de arriba: A
%% Ordeno por X:

[coord2A,indiceXA] = sort(nodosEle(IndexCaraArribaZ,1));

IndexOrdXA = IndexCaraArribaZ(indiceXA);
IndexPrimerosA = IndexOrdXA(1:2);
IndexSegundosA = IndexOrdXA(3:4);

%% Ordeno por Y:
[coord3PA,indiceYP] = sort(nodos(IndexPrimerosA,2));
[coord3SA,indiceYS] = sort(nodos(IndexSegundosA,2));

IndexOrdAP = IndexPrimerosA(indiceYP);
IndexPrimeroAP = IndexOrdAP(2);
IndexSegundoAP = IndexOrdAP(1);

IndexOrdAS = IndexSegundosA(indiceYS);
IndexPrimeroAS = IndexOrdAS(2);
IndexSegundoAS = IndexOrdAS(1);

CaraArribaOrd = [IndexPrimeroAP IndexSegundoAP IndexPrimeroAS IndexSegundoAS]; % Nodos ordenados de la cara superior.

%% Para la cara de abajo: B

%% Ordeno por X:

[coord2B,indiceXB] = sort(nodosEle(IndexCaraAbajoZ,1));

IndexOrdXB = IndexCaraAbajoZ(indiceXB);
IndexPrimerosB = IndexOrdXB(1:2);
IndexSegundosB = IndexOrdXB(3:4);

%% Ordeno por Y:
[coord3PB,indiceYP] = sort(nodosEle(IndexPrimerosB,2));
[coord3SB,indiceYS] = sort(nodosEle(IndexSegundosB,2));

IndexOrdBP = IndexPrimerosB(indiceYP);
IndexPrimeroBP = IndexOrdBP(2);
IndexSegundoBP = IndexOrdBP(1);

IndexOrdBS = IndexSegundosB(indiceYS);
IndexPrimeroBS = IndexOrdBS(1);
IndexSegundoBS = IndexOrdBS(2);

CaraAbajoOrd = [IndexPrimeroBP IndexSegundoBP IndexPrimeroBS IndexSegundoBS]; % Nodos ordenados de la cara inferior.

%% Elementos ordenados:
IndexOrd = [CaraArribaOrd CaraAbajoOrd];
elementosOrd(iele,:) =  elementos(iele,IndexOrd);

% coord = [coord1 coord2A coord3PA coord3SA coord2B coord3PB coord3SB];

end

end

