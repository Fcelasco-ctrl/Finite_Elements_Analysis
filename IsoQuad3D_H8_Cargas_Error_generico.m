tic
clear
clc
close all

%% 1 - Peticiones de salida:

eleType = 'H8FACU'; % Segun como sean deducidas las shapefunds

% Graficos:
scale = 400000; % Escala de magnificación de desplazamientos. 100 para q
show = 0; % 0; Para mostrar numeros de nodos en el meshplot.
  % Colores:
defocolor = 'r';
meshcolor = 'b';
inferiorcolor = 'b';
unioncolor = 'b'; % 'g'
superiorcolor = 'b'; % 'k'
sopcolor = 'b'; % 'm'

% Cargas en la estructura:
q = -0.0095; % Carga distribuida [N/mm^2] % -0.18
Rpunt = 0; % Cargas puntuals en la punta izquierda inferior % -10000 [N] (Test)

% Elegir la cantidad de superelementos a usar:
NsuperElemLineal = 10; % Por default 10
NsuperElemRotados = 12; % Por default 12

% Tension a graficar:

S = 1; 
%% 2 - Discretizacion: 3D f(x,y,z)
 
% Carga de mallas: 

% Malla semilla "alivianamientos" superelem1.txt tiene el (0,0,0) en el centro
% de la malla:

nodos = load('nodos_superelem1.txt'); 
elementos = load('elementos_superelem1.txt');

% Malla semilla "union" union.txt tiene el (0,0,0) en la longitud total de
% la cercha:

nodosUn = load('nodos_union_remesh.txt');
elementosUn = load('elementos_union_remesh.txt');

% Tomo los datos que me sirven:
% Malla alivianamientos:
nodos = nodos(:,2:4);
elementos = elementos(:,2:9);
[elementosOrd] = Elements_sorter(nodos,elementos);

elementos = elementosOrd;
nodosPlot = nodos; % Nodos para ploteo de tensiones.

% Malla union:
nodosUn = nodosUn(:,2:4);
elementosUn = elementosUn(:,2:9);

elementosUnion = elementosUn;

% Malla semilla "soporte" .txt tiene el (0,0,0) en la altura total de
% la cercha:

nodosSop = load('nodos_soporte.txt');
elementosSop = load('elementos_soporte.txt');

% Tomo los datos que me sirven:

nodosSop = nodosSop(:,2:4);
elementosSop = elementosSop(:,2:9);

%% 3 - Propiedades del Material:
E=200e6;
nu=0.3;

% Constitutivo 3D:
A = ((1-nu)*E)/((1+nu)*(1-2*nu));
B = (nu*E)/((1+nu)*(1-2*nu));
G = E/(2*(1+nu));
C = [A B B 0 0 0
     B A B 0 0 0
     B B A 0 0 0
     0 0 0 G 0 0
     0 0 0 0 G 0
     0 0 0 0 0 G];

%% 3a - Definiciones de mallas y garficos de mallas semilla:
%---------------------------------------------
% 3a - Definiciones: Malla semilla "ALIVIANAMIENTOS"
%---------------------------------------------
nDofNod = 3;                    % grados de libertad por nodo
nel = size(elementos,1);         % cantidad de elementos
nNod = size(nodos,1);           % cantidad de nodos
nNodEle = size(elementos,2);     % nodos por elemento
nDofTot = nDofNod*nNod;         % grados de libertad de malla semilla
nDims = size(nodos,2);          % dimensiones del problema
nodeDoscale = reshape(1:nDofTot,nDofNod,nNod)'; % Matriz de grados de libertad
dof = nodeDoscale; % Matriz de dof por nodos

%---------------------------------------------
% 3b - Definiciones: Malla semilla "UNION"
%---------------------------------------------
nDofNodUn = 3;                    % grados de libertad por nodo
nelUn = size(elementosUnion,1);         % cantidad de elementos
nNodUn = size(nodosUn,1);           % cantidad de nodos
nNodEleUn = size(elementosUnion,2);     % nodos por elemento
nDofTotUn = nDofNodUn*nNodUn;         % grados de libertad de malla semilla
nDimsUn = size(nodosUn,2);          % dimensiones del problema
nodeDoscaleUn = reshape(1:nDofTotUn,nDofNodUn,nNodUn)'; % Matriz de grados de libertad
dofUn = nodeDoscaleUn; % Matriz de dof por nodos

%---------------------------------------------
% 3C - Definiciones: Malla semilla "SOPORTE"
%---------------------------------------------

nDofNodSop = 3;                    % grados de libertad por nodo
nelSop = size(elementosSop,1);         % cantidad de elementos
nNodSop = size(nodosSop,1);           % cantidad de nodos
nNodEleSop = size(elementosSop,2);     % nodos por elemento
nDofTotSop = nDofNodSop*nNodSop;         % grados de libertad de malla semilla
nDimsSop = size(nodosSop,2);          % dimensiones del problema
nodeDoscaleSop = reshape(1:nDofTotSop,nDofNodSop,nNodSop)'; % Matriz de grados de libertad
dofSop = nodeDoscaleSop; % Matriz de dof por nodos

% Grafico de malla en 3D:

figure()
Meshplot(elementos,nodos,meshcolor,show) 
title('Malla semilla: Perfil con alivianamientos');

figure()
Meshplot(elementosUn,nodosUn,meshcolor,show) 
title('Malla: Unión');

figure()
Meshplot(elementosSop,nodosSop,meshcolor,show) 
title('Malla: Soporte');

%% 4 - Matrices de rigidez: 
% Puntos de Gauss: Para 3D
a = 1/sqrt(3);
%        ksi eta zeta
upg = a*[-1  -1   1
         -1   1   1
          1  -1   1
          1   1   1
         -1  -1  -1
         -1   1  -1
          1  -1  -1
          1   1  -1];
      
npg = size(upg,1);
wpg = ones(1,8);
%----------------------------------------------------------------
%% 4a - Matriz de rigidez Global: MALLA SEMILLA "ALIVIANAMIENTOS":
%----------------------------------------------------------------
tic
K = zeros(nDofTot);
cont = 0; % Chequeo de entrada a ciclo.
for iele = 1:nel
    Ke = zeros(nDofNod*nNodEle);
    nodosEle = nodos(elementos(iele,:),:);
    for ipg = 1:npg
        
        % Puntos de Gauss:
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
 
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac = dN*nodosEle;
        
        % Derivadas de las funciones de forma respecto de x,y,z
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        B = zeros(size(C,2),nDofNod*nNodEle);
        B(1,1:3:nDofNod*nNodEle-2) = dNxy(1,:);
        B(2,2:3:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,3:3:nDofNod*nNodEle) = dNxy(3,:);
        B(4,1:3:nDofNod*nNodEle-2) = dNxy(2,:);
        B(4,2:3:nDofNod*nNodEle-1) = dNxy(1,:);
        B(5,2:3:nDofNod*nNodEle-1) = dNxy(3,:);
        B(5,3:3:nDofNod*nNodEle) = dNxy(2,:);
        B(6,1:3:nDofNod*nNodEle-2) = dNxy(3,:);
        B(6,3:3:nDofNod*nNodEle) = dNxy(1,:);
        
        Djac = det(jac);
        Ke = Ke + B'*C*B*wpg(ipg)*Djac;
        cont = cont + 1;
    end
    
    eleDoscale = nodeDoscale(elementos(iele,:),:);
    eleDoscale = reshape(eleDoscale',[],1);
    K(eleDoscale,eleDoscale) = K(eleDoscale,eleDoscale) + Ke; 
end
% Chequeo de autovalores:
eigenValK = eig(K);

% Chequeo de ensamble y Chqueo de diagonalidad:
if max(max(abs(K-K')))< 9e-3 && trace(K) > 0 && abs(sum(eigenValK(size(eigenValK,1)-5:size(eigenValK,1)))) < 1e-4
   disp('El ensamble de K es correcto');
   disp('---------------------------------');
%    figure()
%    spy(K);
%    title('Matriz de rigidez KLineal');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KLineal MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');
% rcond(K);

%----------------------------------------------------------------
%% 4b - Matriz de rigidez Global: SUB MALLA "UNION":
%----------------------------------------------------------------
tic
KUn = zeros(nDofTotUn);
cont = 0; % Chequeo de entrada a ciclo.
for iele = 1:nelUn
    KeUn = zeros(nDofNodUn*nNodEleUn);
    nodosEleUn = nodosUn(elementosUnion(iele,:),:);
    for ipg = 1:npg
        
        % Puntos de Gauss:
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
 
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac = dN*nodosEleUn;
        
        % Derivadas de las funciones de forma respecto de x,y,z
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        B = zeros(size(C,2),nDofNodUn*nNodEleUn);
        B(1,1:3:nDofNodUn*nNodEleUn-2) = dNxy(1,:);
        B(2,2:3:nDofNodUn*nNodEleUn-1) = dNxy(2,:);
        B(3,3:3:nDofNodUn*nNodEleUn) = dNxy(3,:);
        B(4,1:3:nDofNodUn*nNodEleUn-2) = dNxy(2,:);
        B(4,2:3:nDofNodUn*nNodEleUn-1) = dNxy(1,:);
        B(5,2:3:nDofNodUn*nNodEleUn-1) = dNxy(3,:);
        B(5,3:3:nDofNodUn*nNodEleUn) = dNxy(2,:);
        B(6,1:3:nDofNodUn*nNodEleUn-2) = dNxy(3,:);
        B(6,3:3:nDofNodUn*nNodEleUn) = dNxy(1,:);
        
        Djac = det(jac);
        KeUn = KeUn + B'*C*B*wpg(ipg)*Djac;
        cont = cont + 1;
    end
    
    eleDoscaleUn = nodeDoscaleUn(elementosUnion(iele,:),:);
    eleDoscaleUn = reshape(eleDoscaleUn',[],1);
    KUn(eleDoscaleUn,eleDoscaleUn) = KUn(eleDoscaleUn,eleDoscaleUn) + KeUn; 
end
% Chequeo de autovalores:
eigenValKUn = eig(KUn);
% Chequeo de ensamble:
if max(max(abs(KUn-KUn')))< 9e-3 && trace(KUn) > 0 && abs(sum(eigenValKUn(size(eigenValKUn,1)-5:size(eigenValKUn,1)))) < 1e-4
   disp('El ensamble de KUn es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KUn);
%    title('Matriz de rigidez KUn');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KUn MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');
% rcond(KUn);

%----------------------------------------------------------------
%% 4c - Matriz de rigidez Global: SUB MALLA "SOPORTE":
%----------------------------------------------------------------
tic
KSop = zeros(nDofTotSop);
cont = 0; % Chequeo de entrada a ciclo.
for iele = 1:nelSop
    KeSop = zeros(nDofNodSop*nNodEleSop);
    nodosEleSop = nodosSop(elementosSop(iele,:),:);
    for ipg = 1:npg
        
        % Puntos de Gauss:
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
 
        % Derivadas de las funciones de forma respecto de ksi, eta
        dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac = dN*nodosEleSop;
        
        % Derivadas de las funciones de forma respecto de x,y,z
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        B = zeros(size(C,2),nDofNodSop*nNodEleSop);
        B(1,1:3:nDofNodSop*nNodEleSop-2) = dNxy(1,:);
        B(2,2:3:nDofNodSop*nNodEleSop-1) = dNxy(2,:);
        B(3,3:3:nDofNodSop*nNodEleSop) = dNxy(3,:);
        B(4,1:3:nDofNodSop*nNodEleSop-2) = dNxy(2,:);
        B(4,2:3:nDofNodSop*nNodEleSop-1) = dNxy(1,:);
        B(5,2:3:nDofNodSop*nNodEleSop-1) = dNxy(3,:);
        B(5,3:3:nDofNodSop*nNodEleSop) = dNxy(2,:);
        B(6,1:3:nDofNodSop*nNodEleSop-2) = dNxy(3,:);
        B(6,3:3:nDofNodSop*nNodEleSop) = dNxy(1,:);
        
        Djac = det(jac);
        KeSop = KeSop + B'*C*B*wpg(ipg)*Djac;
        cont = cont + 1;
    end
    
    eleDoscaleSop = nodeDoscaleSop(elementosSop(iele,:),:);
    eleDoscaleSop = reshape(eleDoscaleSop',[],1);
    KSop(eleDoscaleSop,eleDoscaleSop) = KSop(eleDoscaleSop,eleDoscaleSop) + KeSop; 
end
% Chequeo de autovalores:
eigenValKSop = eig(KSop);
% Chequeo de ensamble:
if max(max(abs(KSop-KSop')))< 9e-3 && trace(KSop) > 0 && abs(sum(eigenValKSop(size(eigenValKSop,1)-5:size(eigenValKSop,1)))) < 1e-4
   disp('El ensamble de KSop es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KSop);
%    title('Matriz de rigidez KSop');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KSop MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');
% rcond(KSop);


%% 5 - Condensación del Super Elementos ALIVIANAMIENTOS Lineales y Rotados, y ordenamiento KUn y KSop

% Nota 1 : SE DISTINGUE ENTRE las sub mallas de ALIVIANAMIENTOS "LINEALES" Y LOS "ROTADOS":
% Nota 2 : La matriz del super elemento "propagada": Nint veces. corresponde a 
% el super elemento que mide: 200 [mm] de largo

% Ordeno los nodosB de la derecha y de la izquierda segun un mismo 
% criterio, para poder pegar bien las mallas del super elemento propagadas:
% Se usa el criterio de "abajo para arriba" y "desde el fondo hacia
% adelante".

% Matriz de rigidez del super elemento "alivianamientos" semilla: Condensación. Nodos externos (b) e internos (i):
% Nodos externos (b),son los que estan en los laterales de la malla semilla respecto del centro
% x = 0 y x = -100, Los dof de la Kbmm (SE) se ordenan como: dofI(b)
% dofD(b). En la descondensación se encuentran los dof(i).

%--------------------------------------------
%% 5a - Para la sub malla semilla ALIVIANAMIENTOS:
%--------------------------------------------
tic
%--------------------------
% 5a - Nodos Boundary (b): 
%--------------------------

% 5a- 1 Cuantos hay en las puntas, techo y piso del perfil  exteriores:
cb = 0;
cbi = 0; % Nodos (b) izquierdos externos
cbd = 0; % Nodos (b) derechos externos
cbt = 0; % Nodos (b) del techo externos: No lo uso
for inod = 1:nNod
            if abs(nodos(inod,1)+100)<1e-9 
               cb = cb +1;
               cbi = cbi + 1;
            else
                if abs(nodos(inod,1)-100)<1e-9 
                 cb = cb +1;
                 cbd = cbd + 1;
                end
            end
end

% 5a- 2 Busco los nodos del techo (ala) de mi malla:
% Cuantos hay?:
for inod = 1:nNod
 if abs(nodos(inod,3)-100)<1e-9
    cbt = cbt + 1;
 end
end

% 5a- 3 Creo el vector nodosB, con los indices de los nodos externos:
nodosB = zeros(cb,1);
nodosBi = zeros(cbi,1); % Indice de nodos (b) izquierdos externos.
nodosBd = zeros(cbd,1); % Indice de nodos (b) derechos externos.
posB = 1;
posBi = 1;
posBd = 1;
for inod = 1:nNod
            if abs(nodos(inod,1)-100)<1e-9
               nodosB(posB) = inod;
               nodosBi(posBi) = inod;
               posB = posB + 1;
               posBi = posBi + 1;
            else
                if abs(nodos(inod,1)+100)<1e-9
                   nodosB(posB) = inod;
                   nodosBd(posBd) = inod;
                   posB = posB + 1;
                   posBd = posBd + 1;
                end
            end
end 

% 5a - 4 Creo el vector nodosBt, con los indices de los nodos externos del techo: 
% Busco los nodos de los elementos de la cara de arriba: z = 100.   

% nodosBt = zeros(cbt,1); % Nodos del techo
% posBt = 1;
% for inod = 1:nNod
%     if abs(nodos(inod,3)-100)<1e-9 
%         nodosBt(posBt) = inod;
%         posBt = posBt + 1;
%     end
% end

%-------------------------
% 5b - Nodos Interior (i):
%-------------------------

%Todos los que no son nodosB.

% 5b - 1 ¿Cuantos hay?:

ci = size(nodos,1) - size(nodosB,1);

% 5b - 2 Creo el vector nodosI, con los indices de los nodos internos:
posI = true(nNod,1);
posI(nodosB,1) = false;

nodosI = find(posI == 1); % Nodos internos chequados: Los que no se conectan.

%----------------------------
% 5c - Ordenamiento de nodos:
%----------------------------

% Ordeno los nodosB de la derecha y de la izquierda segun un mismo 
% criterio, para poder pegar bien las mallas del super elemento propagadas:

% 5c - 1 Ordeno la cara izquierda cI:
ordenZ = unique(sort(nodos(nodosBi,3)));
cIo =zeros(size(ordenZ,1),length(nodosBi));
for i=1:length(ordenZ)
  zz = find(nodos(nodosBi,3)==ordenZ(i)); %nodos por altura en Z
  ordenY = unique(sort(nodos(nodosBi(zz),2)));
  for j=1:length(zz)
   yy(j) = find(nodos(nodosBi(zz),2)==ordenY(j));
  end
   cIo(i,[1:length(zz(yy))])=nodosBi(zz(yy));
   yy = 0;
end

cIor=reshape(cIo',[],1); % Indices de nodosI ordenados.

% Eliminando filas y columnas de ceros:
BC=sum(abs(cIor),1)==0;  
cIor(sum(abs(cIor),2)==0,:)=[ ];
cIor(:,BC)=[ ];
nodosBi=cIor;

% 5c - 2 Ordeno la cara derecha cD:
ordenZ = unique(sort(nodos(nodosBd,3)));
cDo =zeros(size(ordenZ,1),length(nodosBd));
for i=1:length(ordenZ)
  zz = find(nodos(nodosBd,3)==ordenZ(i)); %nodos por altura en Z
  ordenY = unique(sort(nodos(nodosBd(zz),2)));
  for j=1:length(zz)
   yy(j) = find(nodos(nodosBd(zz),2)==ordenY(j));
  end
   cDo(i,[1:length(zz(yy))])=nodosBd(zz(yy));
   yy = 0;
end

cDor=reshape(cDo',[],1); % Indices de nodosD ordenados.

% Elimino filas y columnas de ceros
BC=sum(abs(cDor),1)==0;  
cDor(sum(abs(cDor),2)==0,:)=[ ];
cDor(:,BC)=[ ];
nodosBd=cDor;

% 5c - 3 Renombrando variables para la parte de la matriz de rigidez:
nodosBii = nodosBi; % Nodos izquierdos, para buscar los dofI(b)
nodosBdd = nodosBd; % Nodos derechos, para buscar los dofD(b)

%--------------------------------
%% 5b - Para la sub malla  UNION:
%--------------------------------
% Ordeno los nodos que se conectan con la viga de abajo que se unirian 
% con nodosBdd del final de la matriz de rigidez propagada de
% ALIVIANAMIENTOS Lineales, con el l mismo criterio, para poder pegar 
% bien las mallas del super elemento propagadas:

%--------------------------
% 5a - Nodos Boundary (b): 
%--------------------------

% 5a - 1 Cuantos nodos hay en las uniones de la UNION:
% y en los nodos que toman BC:
cbUn = 0;
cbUnBC = 0; % Nodos (b) que toman boundary conditions del ensamble final.
cbUinf = 0; % Nodos (b) izquierdos conectores con viga inferior-
cbUsup = 0; % Nodos (b) izquierdos conectores de la viga superior
            % z = mx + b Recta que contiene a nodos conectores viga superior
            % Esta contenidos si z - mx - b = 0
m = -1.4295; % Si no cambia la geometría de Malla Semilla UNION
b = -2.5553e+03; % Si no cambia la geometría de Malla Semilla UNION
for inod = 1:nNodUn
            if abs(nodosUn(inod,1)+2000)<1e-9 && abs(nodosUn(inod,3))<=100 
              cbUn = cbUn+ 1;
               cbUinf = cbUinf + 1;
            else
                if abs(nodosUn(inod,3)-m*nodosUn(inod,1)-b)<0.5
                 cbUn = cbUn + 1;
                 cbUsup = cbUsup + 1;
            else
                if abs(nodosUn(inod,1)+2200)<1e-9 && abs(nodosUn(inod,3))<=100
                cbUnBC = cbUnBC + 1;  
                    
                end
                end
            end
end

% 5a- 2 Creo el vector nodosBUn, con los indices de los nodos externos de la UNION:
nodosUnTot = 1:1:max(max(elementosUn)); % Indice de nodos totales de la malla UNION
nodosBUn = zeros(cbUn,1);
nodosBinfUn = zeros(cbUinf,1); % Indice de nodos (b) de conexion viga inferior externos.
nodosBsupUn = zeros(cbUsup,1); % Indice de nodos (b) de conexion de viga superior externos.
nodosBUnBC = zeros(cbUnBC,1); % Indice de nodos (b) que toman condiciones de borde de todo el modelo.
posBUn = 1;
posBinfUn = 1;
posBsupUn = 1;
posBUnBC = 1; 
for inod = 1:nNodUn
            if abs(nodosUn(inod,1)+2000)<1e-9 && abs(nodosUn(inod,3))<=100 
               nodosBUn(posBUn) = inod;
               nodosBinfUn(posBinfUn) = inod;
               posBUn = posBUn + 1;
               posBinfUn = posBinfUn + 1;
            else
                if abs(nodosUn(inod,3)-m*nodosUn(inod,1)-b)<0.5
                   nodosBUn(posBUn) = inod;
                   nodosBsupUn(posBsupUn) = inod;
                   posBUn = posBUn + 1;
                   posBsupUn = posBsupUn + 1;
            else
                if abs(nodosUn(inod,1)+2200)<1e-9 && abs(nodosUn(inod,3))<=100
                   nodosBUnBC(posBUnBC) = inod; 
                   posBUnBC = posBUnBC + 1;
                end
                end
            end
end

%----------------------------
% 5b - Ordenamiento de nodos:
%----------------------------

% Ordeno los nodosB de la derecha y de la izquierda segun un mismo 
% criterio, para poder pegar bien las mallas del super elemento propagadas:

% 5b - 1 Ordeno los de la conexion viga INFERIOR externos: Criterio
% IZQUIERDA
ordenZUn = unique(sort(nodosUn(nodosBinfUn,3)));
cIoUn =zeros(size(ordenZUn,1),length(nodosBinfUn));
for i=1:length(ordenZUn)
  zz = find(nodosUn(nodosBinfUn,3)==ordenZUn(i)); %nodos por altura en Z
  ordenY = unique(sort(nodosUn(nodosBinfUn(zz),2)));
  for j=1:length(zz)
   yy(j) = find(nodosUn(nodosBinfUn(zz),2)==ordenY(j));
  end
   cIoUn(i,[1:length(zz(yy))])=nodosBinfUn(zz(yy));
   yy = 0;
end

cIorUn=reshape(cIoUn',[],1); % Indices de nodos inferiores externos ordenados.

% Eliminando filas y columnas de ceros:
BC=sum(abs(cIorUn),1)==0;  
cIorUn(sum(abs(cIorUn),2)==0,:)=[ ];
cIorUn(:,BC)=[ ];
nodosBinfUn=cIorUn;

% 5b - 2 Ordeno los de la conexion viga SUPERIOR externos: Criterio
% IZQUIERDA
ordenZUnSup = unique(sort(nodosUn(nodosBsupUn,3)));
cIoUnSup =zeros(size(ordenZUnSup,1),length(nodosBsupUn));
for i=1:length(ordenZUnSup)
  zz = find(nodosUn(nodosBsupUn,3)==ordenZUnSup(i)); %nodos por altura en Z
  ordenY = unique(sort(nodosUn(nodosBsupUn(zz),2)));
  for j=1:length(zz)
   yy(j) = find(nodosUn(nodosBsupUn(zz),2)==ordenY(j));
  end
   cIoUnSup (i,[1:length(zz(yy))])= nodosBsupUn(zz(yy));
   yy = 0;
end

cIorUnSup=reshape(cIoUnSup',[],1); % Indices de nodos inferiores externos ordenados.

% Eliminando filas y columnas de ceros:
BC=sum(abs(cIorUnSup),1)==0;  
cIorUnSup(sum(abs(cIorUnSup),2)==0,:)=[ ];
cIorUnSup(:,BC)=[ ];
nodosBsupUn=cIorUnSup;

% 5c - 2 Renombrando variables para la parte de la matriz de rigidez:
nodosBUnSup = nodosBsupUn; % Nodos superiores  de la union ordenados criterio IZQUIERDA
nodosBUnInf = nodosBinfUn; % Nodos inferiores  de la union ordenados criterio IZQUIERDA
nodosBUnbc = nodosBUnBC;   % Nodos inferiores (b) que reciben bc(algo)
nodNoInt = [nodosBUnSup nodosBUnInf nodosBUnbc];
nodosUnInt = setdiff(nodosUnTot,nodNoInt)';

%----------------------------------
%% 5c - Para la sub malla  SOPORTE:
%----------------------------------

%--------------------------
% 5a - Nodos Boundary (b): 
%--------------------------

% 5a - 1 Cuantos nodos hay en las uniones de SOPORTE:
% y en los nodos que toman BC:
cbSopBC = 0;  % Nodos (b) que toman boundary conditions del ensamble final.
cbSopCon = 0; % Nodos (b) derechos conectores de la viga superior
              % z = mx + b Recta que contiene a nodos conectores viga superior
              % Esta contenidos si z - mx - b = 0
mSop = -1.4279;  % Si no cambia la geometría de Malla Semilla UNION
bSop = 1631.7;   % Si no cambia la geometría de Malla Semilla UNION
for inod = 1:nNodSop
    
            if abs(nodosSop(inod,1)+12.8)<1e-9 
              cbSopBC = cbSopBC + 1;
            end
            
            if abs(nodosSop(inod,3)-mSop*nodosSop(inod,1)-bSop)<0.5
                 cbSopCon = cbSopCon + 1;
            end
          
end

nodosSopTot = 1:1:max(max(elementosSop)); % Indice de nodos totales de la malla SOPORTE
nodosBSopCon = zeros(cbSopCon,1); % Indice de nodos (b) de conexion viga superior externos.
nodosBSopBC = zeros(cbSopBC,1); % Indice de nodos (b) que toman condiciones de borde de todo el modelo.
posBSop = 1;
posBSopCon = 1;
posBSopBC = 1; 
for inod = 1:nNodSop
            if abs(nodosSop(inod,1)+12.8)<1e-9  
               nodosBSopBC(posBSopBC) = inod;
               posBSopBC = posBSopBC + 1;
            end
            
            if abs(nodosSop(inod,3)-mSop*nodosSop(inod,1)-bSop)<0.5
               nodosBSopCon(posBSopCon) = inod;
               posBSopCon = posBSopCon + 1;
            end     
end

% 5b - 2 Ordeno los de la conexion viga SUPERIOR externos: Criterio
% IZQUIERDA
ordenZSop = unique(sort(nodosSop(nodosBSopCon,3)));
cIoSop = zeros(size(ordenZSop,1),length(nodosBSopCon));
for i=1:length(ordenZSop)
  zz = find(nodosSop(nodosBSopCon,3)==ordenZSop(i)); %nodos por altura en Z
  ordenY = unique(sort(nodosSop(nodosBSopCon(zz),2)));
  for j=1:length(zz)
   yy(j) = find(nodosSop(nodosBSopCon(zz),2)==ordenY(j));
  end
   cIoSop(i,[1:length(zz(yy))])= nodosBSopCon(zz(yy));
   yy = 0;
end

cIorSop=reshape(cIoSop',[],1); % Indices de nodos inferiores externos ordenados.

% Eliminando filas y columnas de ceros:
BC=sum(abs(cIorSop),1)==0;  
cIorSop(sum(abs(cIorSop),2)==0,:)=[ ];
cIorSop(:,BC)=[ ];

% 5c - 3 Renombrando variables para la parte de la matriz de rigidez:
nodosBSopConOrd = cIorSop;
nodosBSopbc = nodosBSopBC;
nodosSopInterv = [nodosBSopConOrd' nodosBSopbc'];
nodosSopInt = setdiff(nodosSopTot,nodosSopInterv)';

%--------------------------------------
%% 5d - "Explicit matrix operations": Kbb chapter 10 FEA paper.
%--------------------------------------

% Matriz de rigidez Kbb condensada by "Explicit matrix operations": 
% Super Elemento "Alivianamientos" de una malla semilla cargada.
% Matriz de rigidez local del Super Elemento ya ordenada simetricamente 
% con sus dofI y dofD con el mismo orden y criterio de selección:

% 5d - 1 Reshape de los dofs:
boundDofs = dof(nodosI,:); % Dofs biundary sin hacer el reshape: U V W
boundaryDofs = reshape(dof([nodosBii nodosBdd],:)',[],1); % dof de los nodos externos
interiorDofs = reshape(dof(nodosI,:)',[],1); % dof de los nodos internos para una malla semilla
boundaryDofsI = reshape(dof(nodosBii,:)',[],1); % dof boundary izq ORDENADOS
boundaryDofsD = reshape(dof(nodosBdd,:)',[],1); % dof boundatry derechos ORDENADOS
nDofI = size(boundaryDofsI,1); % numero de dof cara izq. (Es indistinto, se repite la misma cantidad en ambas caras)
nDofD = size(boundaryDofsD,1); % numero de dof cada derecha
nDofBoundaryFace = nDofI; % dof totales en una cara del super elemento. (Es igual para ambas caras)
nDofBoundary = nDofI + nDofD; % dof totales en los extremos del super elemento.

% 5d - 2 Matriz de rigidez del super elemento "ALIVIANAMIENTOS Lineales":

% Matriz de rigidez ordenada como ub vb wb, segun: nodosBii nodosBdd.
Kbbm = K(boundaryDofs,boundaryDofs) - (K(boundaryDofs,interiorDofs)*inv(K(interiorDofs,interiorDofs))*K(interiorDofs,boundaryDofs)); % Matriz  de rigidez del super elemento FACU
disp('Tiempo de condensado');
toc
% 5d - 3 Ensamble de super elementos:

%------------------------------------------------------------
%% 5e - Sub Malla ALIVIANAMIENTOS Lineales: Izq + (Int) + Der
%------------------------------------------------------------
tic
% NsuperElemLineal = 10; % Cantidad de super elementos a ensamblar

nint = NsuperElemLineal + 1; % Numero de interfaces de nodos resultantes del ensamble

% Inicializo matriz de ceros para ensamblar:
KpropLineal = zeros(nint*nDofBoundaryFace,nint*nDofBoundaryFace);

% Ensamblo:
for idof = 1:NsuperElemLineal
    dofSEa = nDofBoundaryFace*(idof-1) + 1;
    dofSEb = nDofBoundaryFace*(idof+1);
    dofSElin = dofSEa:dofSEb;
    KpropLineal(dofSElin,dofSElin) = KpropLineal(dofSElin,dofSElin) + Kbbm;
end
% Chequeo de ensamble:
if max(max(abs(KpropLineal-KpropLineal')))< 9e-3 && trace(KpropLineal) > 0
   disp('El ensamble de KpropLineal es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KpropLineal);
%    title('K de sub malla de Super Elementos Lineales Propagada');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KpropLineal MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

nDofBoundLineal = size(KpropLineal,1); % Numero de dof boundary TOTALES de mi malla de super elementos.
nDofInteriorLineal = size(interiorDofs,1); % Numero de dof interiores de UNA malla semilla.

%-----------------------------------------
%% 5f - Sub Malla ALIVIANAMIENTOS Rotados:
%-----------------------------------------
tic
% NsuperElemRotados = 12; % Cantidad de super elementos a ensamblar

nintRot = NsuperElemRotados + 1; % Numero de interfaces de nodos resultantes del ensamble

% Rotación Pytch: Armo matriz T del elemento: Bibliografia "Directional
% Cosine Matrices". pag 1.
theta = 35;
Lamda = [ cosd(theta) 0  -sind(theta)
              0       1       0
          sind(theta) 0  cosd(theta)];

T = zeros(size(K));

cont = [1 2 3];
sumador = 3*ones(1,3);
for Nod = 1:nNod
    T(cont,cont) = Lamda; 
    cont = cont + sumador;
end

Krot = T'*K*T; % Matriz de rigidezde mi malla semilla "ALIVIANAMIENTOS" rotada.
if max(max(abs(Krot-Krot')))< 9e-3 && trace(Krot) > 0
   disp('La obtención de Krot es correcta');
   disp('---------------------------------');
%    figure()
%    spy(Krot);
%    title('K de sub malla de Super Elementos Rotados');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ Krot MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

eigenKrot = eig(Krot);

% Condensación para "Alivianamientos rotados":
KbbmRot = Krot(boundaryDofs,boundaryDofs) - (Krot(boundaryDofs,interiorDofs)*inv(Krot(interiorDofs,interiorDofs))*Krot(interiorDofs,boundaryDofs)); 

% Inicializo matriz de ceros para ensamblar:
KpropRotados = zeros(nintRot*nDofBoundaryFace,nintRot*nDofBoundaryFace);

% Ensamblo y roto:
for idof = 1:NsuperElemRotados
    dofSEa = nDofBoundaryFace*(idof-1) + 1;
    dofSEb = nDofBoundaryFace*(idof+1);
    dofSE = dofSEa:dofSEb;
%     KbbmRot = T'*Kbbm*T;
    KpropRotados(dofSE,dofSE) = KpropRotados(dofSE,dofSE) + KbbmRot;
end
% Chequeo de ensamble:
if max(max(abs(KpropRotados-KpropRotados')))< 9e-3 && trace(KpropRotados) > 0 && abs(sum(eigenKrot(size(eigenKrot,1)-5:size(eigenKrot,1)))) < 1e-4
   disp('El ensamble de KpropRotados es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KpropRotados);
%    title('K de malla de Super Elementos Rotados Propagada');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KpropRotados MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

nDofBoundRotados= size(KpropRotados,1); % Numero de dof boundary TOTALES de mi malla de super elementos.

%----------------------------------------
%% 5g - Sub Malla UNION ORDENADA Inf - bc - Sup:
%----------------------------------------
tic
InferiorUnDofs = reshape(dofUn(nodosBUnInf,:)',[],1)'; % dof de los nodos inferiores ORDENADOS
SuperiorUnDofs = reshape(dofUn(nodosBUnSup,:)',[],1)'; % dof boundary izq ORDENADOS
InteriorDofs = reshape(dofUn(nodosUnInt,:)',[],1)'; % dof de los nodos inferiores ORDENADOS
InteriorDofsBC = reshape(dofUn(nodosBUnbc,:)',[],1)'; % dof de los nodos inferiores ORDENADOS

ordenUn = [InferiorUnDofs InteriorDofs InteriorDofsBC SuperiorUnDofs];
KUnOrd = KUn(ordenUn,ordenUn);
eigenKUnOrd = eig(KUnOrd);

% Cantidad de dof con las denominaciones anteriores: 
nDofInfUn = size(InferiorUnDofs,2);
nDofIntUn = size(InteriorDofs,2);
nDofIntBCUn = size(InteriorDofsBC,2);
nDofSupUn = size(SuperiorUnDofs,2);
nDofTotUn = nDofInfUn + nDofIntUn + nDofIntBCUn + nDofSupUn; % Nodos totales

% Chequeo de ensamble:
if max(max(abs(KUnOrd-KUnOrd')))< 9e-3 && trace(KUnOrd) > 0 && abs(sum(eigenKUnOrd(size(eigenKUnOrd,1)-5:size(eigenKUnOrd,1)))) < 1e-4
   disp('El ensamble de KUnOrd es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KUnOrd)
%    title('K sub malla Kun Ordenada');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KUnOrd MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

%--------------------------------------------------
%% 5h - Sub Malla SOPORTE ORDENADA: Con + Int + BC
%--------------------------------------------------
tic
nodosConectoresSop = nodosBSopConOrd; % Nodos conectores a la estructura De la viga SUPERIOR
% nodosConBCSop = nodosBSopbc; % Nodos que toman condiciones de borde
nodosConBCSop = [9 10 11 12 17 18]; % Nodos que toman condiciones de borde
nodosInternosSop = nodosSopInt; % Nodos internos de la malla soporte

conectorSopDofsOrd = reshape(dofSop(nodosConectoresSop,:)',[],1)'; % dof que conectan con la viga superior Ordenados criterio IZQUIERDA
boundaryCondSopDofs = reshape(dofSop(nodosConBCSop,:)',[],1)'; % dof que toman BC
interiorSopDofs = reshape(dofSop(nodosInternosSop,:)',[],1)'; % dof de los nodos interiores

ordenSop = [conectorSopDofsOrd interiorSopDofs boundaryCondSopDofs];
KSopOrd = KSop(ordenSop,ordenSop);
eigenKSopOrd = eig(KSopOrd);

% Chequeo de ensamble:
if max(max(abs(KSopOrd-KSopOrd')))< 9e-3 && trace(KSopOrd) > 0 && abs(sum(eigenKSopOrd(size(eigenKSopOrd,1)-5:size(eigenKSopOrd,1)))) < 1e-4
   disp('El ensamble de KSopOrd es correcto');
   disp('---------------------------------');
%    figure()
%    spy(KSopOrd)
%    title('K sub malla KSopOrd Ordenada');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ KSopOrd MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

%% 6 - Ensamble de matrices: KpropLineal + KUn + KpropRotados + KSopOrd. Respetando el ordenamiento de los nodos.
tic
% Emsamble o concatenamiento de todas las matrices de rigidez de las
% partes del modelo que fueron discretizadas previamente:

% 6a - Ensamble de matrices de rigidez:
% Dof totales de cada matriz de rigidez individual:

nDofKpropLineal = size(KpropLineal,1); % Ensamble de SE
nDofKpropRotados = size(KpropRotados,1); % Ensamble de SE
nDofKUn = size(KUnOrd,1); % Ensamble de K comun
nDofKSop = size(KSopOrd,1); % Ensamble de K comun

% Dof totales de cada conección entre matrices de rigidez:

nDofCon = 3*size(nodosBii,1); % Cantidad de dofs interviniendo en las conecciones.
nCon = 3; % Cantidad de conecciones de mallas semilla que hago

% Dof totales que tiene mi matriz de rigidez TOTAL ensamblada:

nDofKtotal = - nCon*nDofCon + nDofKpropLineal + nDofKUn + nDofKpropRotados + nDofKSop; % Suma total de dofs del problema.

% Parámetros:

dofKtotal = reshape(1:nDofKtotal,3,[])';

% Ensamblo:
Ktotal = zeros(nDofKtotal);

% Lista de indices de matrices individuales:
dofEnsUnLin = 1:nDofCon; % Indice de primeros dofs de KUnOrdenada para ensamblar
dofUnOrdResto = 1:nDofKUn; % Indice de dofs KUnOrdenada, despues del ensamle con KpropLineal

% Lista de indices para Ktotal:
dofEnsLinKtotIndex = dofSElin(1,(size(dofSElin,2)/2)+1:end); % Indice de ultimos dofs de KpropLin a ensamnblar con KUnOrd.
dofUnOrd = nDofKpropLineal+1-nDofCon:(nDofKpropLineal+1)+size(KUnOrd,1)-1-nDofCon; % Indice de segundos dof de KUnOrdenada para ensamblar
dofUnOrdRot = max(dofUnOrd)-nDofCon+1:(max(dofUnOrd)-nDofCon)+nDofKpropRotados; % Indice de dofs de Ktotal para ensamblar KpropRotados 
dofRotSop = (max(dofUnOrdRot)+1)-nDofCon:((max(dofUnOrdRot)+1)-nDofCon)+nDofKSop-1; % Indice de dofs de Ktotal apra ensamlar KSopOrd.

% A - Ensamblo KpropLineal:
Ktotal(1:nDofKpropLineal,1:nDofKpropLineal) = Ktotal(1:nDofKpropLineal,1:nDofKpropLineal) + KpropLineal(1:nDofKpropLineal,1:nDofKpropLineal);
% B - Ensamblo KpropLineal + KUnOrd:
Ktotal(dofEnsLinKtotIndex,dofEnsLinKtotIndex) = Ktotal(dofEnsLinKtotIndex,dofEnsLinKtotIndex) + KUnOrd(dofEnsUnLin,dofEnsUnLin);
% C - Ensamblo KUnOrd:
Ksacar = zeros(size(KUnOrd));
Ksacar(dofEnsUnLin,dofEnsUnLin) = Ksacar(dofEnsUnLin,dofEnsUnLin) + KUnOrd(dofEnsUnLin,dofEnsUnLin);
KUnOrdResto = KUnOrd - Ksacar;
Ktotal(dofUnOrd,dofUnOrd) = Ktotal(dofUnOrd,dofUnOrd) + KUnOrdResto(dofUnOrdResto,dofUnOrdResto);
% D - Ensamblo KUnOrd + KpropRot:
Ktotal(dofUnOrdRot,dofUnOrdRot) = Ktotal(dofUnOrdRot,dofUnOrdRot) +  KpropRotados;
% E - Ensamblo KpropRot + KSopOrd:
Ktotal(dofRotSop,dofRotSop) = Ktotal(dofRotSop,dofRotSop) + KSopOrd;
eigenKtotal = eig(Ktotal);

% Chequeo de ensamble y Chequeo de proporcionalidad directa de matriz de
% rigidez total:
if max(max(abs(Ktotal-Ktotal')))< 9e-3 && trace(Ktotal) > 0
   disp('El ensamble de Ktotal es correcto');
   disp('---------------------------------');
%    figure()
%    spy(Ktotal);
%    title('Matriz de rigidez KpropLineal + KUn + KpropRotados + KSopOrd : Modelo completo');
else
   disp('%%%%%%%%%%%%%%%%%% MATRIZ Ktotal MAL ENSAMBLADA %%%%%%%%%%%%%%');
end
toc
disp('     -------------------------     ');

% Parametros de Ktotal:
nDofTotal = size(Ktotal,1);
% Matriz esparsa de Ktotal:
KtotalSparse = sparse(Ktotal);

%%%%SI NO HAY ERRORES DE ENSAMBLE, CONTINUAR%%%%

%------------------------------------------------------------------
% NOTA: Si las proximas lineas no funcionan:
% Chequear mapeo de dofs, aplicacion de cargas y BC.
% Probar con una carga puntual en algun lado y chequear con ADINA.
% Seguis con post pro.
%------------------------------------------------------------------

%% 7 - Cargas: Cargo solo los dof de "ALIVIANAMIENTOS" rotados.

% 7a - Carga distrubuida en el techo para la malla semilla
% "ALIVIANAMIENTOS" (Completa):

%-----------------------------------------------------------------------------
% NOTA:
% Cargo una malla semilla de super elementos "Alivianamientos Lineales", y la
% oriento hacia abajo en la dirección de ZZ global.
% Utilizo la condensación de carga R(ielem,3) entendiendo que esta en la
% misma dirección Z que para el caso de una viga empotrada en voladizo, con
% carga distribuida: No roto la carga.
%-----------------------------------------------------------------------------

Rs = zeros(nNod,nDofNod); % Vector (matriz) de cargas de superficie por dof para "ALIVIANAMIENTOS rotados"

% Elijo los nodos que estan arriba:
% Busco los nodos del techo (ala) de mi malla:
% Cuantos hay?
cbt = 0;
for inod = 1:nNod
 if abs(nodos(inod,3)-100)<1e-9
    cbt = cbt + 1;
 end
end

% Creo el vector nodosBt, con los indices de los nodos externos del techo: Busco los nodos de los elementos de la cara de arriba:
% z = 100.   
nodosBt = zeros(cbt,1); % Nodos del techo
posBt = 1;
for inod = 1:nNod
    if abs(nodos(inod,3)-100)<1e-9 
        nodosBt(posBt) = inod;
        posBt = posBt + 1;
    end
end
nArriba = nodosBt;
elemsArriba = zeros(size(elementos,1),1);
    % Busco los elementos de arriba:        
        for s = 1:length(nArriba)
            elemsArriba = elemsArriba + (sum(~(elementos-nArriba(s))'))';
        end
        elemsArriba = find(elemsArriba);

% Integro en zeta = 1 Techo de mi problema.       
        a = 1/sqrt(3);
%        ksi eta zeta
upg =   [-a  -a   1
         -a   a   1
          a  -a   1
          a   a   1];
        npg = size(upg,1);
        areas = zeros(nel,1); % Para chequear areas que tienen que sumar 1
        vol = 0; % Chequeo volumenes nuevamente.
%  Integro cargas:       
        for iele = 1:nel
            area = 0;
            if sum(iele == elemsArriba) == 1
                nodosEle = nodos(elementos(iele,:),:);
                fsup = zeros(8,1);
                for ipg = 1:npg
                    % Punto de en los nodos
                    ksi = upg(ipg,1);
                    eta = upg(ipg,2);
                    zeta = upg(ipg,3);
                    % Derivadas de x,y,z respecto de ksi, eta, zeta
                    dN = shapefunsder_3D([ksi eta zeta],eleType);
                    Nf = shapefuns_3D([ksi eta zeta],eleType);
                    jac = dN*nodosEle;
                    % Elijo los jacobianos que van segun: Cook pag 229 -
                    % (6.9 - 8): Integro ksi y eta.
                    jacArea = zeros(2,2);
                    jacArea(1,1) = jac(1,1);
                    jacArea(1,2) = jac(1,2);
                    jacArea(2,1) = jac(2,1);
                    jacArea(2,2) = jac(2,2);
                    
                    Q = q*[1;1;1;1;0;0;0;0];
                    
                    fsup = fsup + Nf'*Nf*Q*wpg(ipg)*det(jacArea);
                    area = area+det(jacArea);
                    Djac = det(jac);
                    vol = vol +  Djac*wpg(ipg);
                    
                end
                % Cargo en dirección Z.
                Rs(elementos(iele,:),3) = Rs(elementos(iele,:),3) + fsup;
            end
            areas(iele) = area;
            iele;
        end
AreaTot = sum(areas); % Area total del techo de una malla semilla "ALIVIANAMIENTOS"
CargaTot = abs(q)*AreaTot;
CargaTotSuma = abs(sum(Rs(:,3)));
% Chequeo de cargas:
if CargaTotSuma == CargaTot
   disp(' Carga de superficie bien integrada');
else
   disp('%%%%CARGA DE SUPERFICIE MAL INTEGRADA%%%%'); 
end
% Hago un reshape para que tenga la misma forma que R:        
RsReshaped = reshape(Rs',[],1);

% Distribuidas en la viga SUPERIOR:

R = RsReshaped;

% 7c - Cargas condensadas para un super elemento: El super elemento queda cargado de esta forma:
% Cargas condensadas para un super elemento.
RcondSE = R(boundaryDofs) - (Krot(boundaryDofs,interiorDofs)*inv(Krot(interiorDofs,interiorDofs)))*R(interiorDofs); 

% 7d - Cargas ensambladas para la malla rotada de mallas semilla:
Rtotal = zeros(nintRot*nDofBoundaryFace,1); % Cargas del ensamble de super elementos.
 for idof = 1:NsuperElemRotados
     
    dofSEa = nDofBoundaryFace*(idof-1) + 1;
    dofSEb = nDofBoundaryFace*(idof+1);
    dofSE = dofSEa:dofSEb;
    Rtotal(dofSE) = Rtotal(dofSE) + RcondSE;
 
 end
 

%% 8 - Condiciones de borde: Aca tengo que hacer el mapeo de los dofs que hice en 6:

% 8a - Mapeo de dofs de la matriz de rididez TOTAL ensamblada:
dofsKtotalKpropLineal = 1:size(KpropLineal,1); % Dofs de KlinealProp mapeados en dofs TOTALES
dofsKtotalKUnOrd = (size(KpropLineal,1)-nDofCon)+1:(size(KpropLineal,1)+1-nDofCon)+size(KUnOrd,1)-1;  % Dofs de KUnOrd mapeados en dofs TOTALES
dofsKtotalKpropRotados = ((size(KpropLineal,1)+1-nDofCon)+size(KUnOrd,1)-1)+1-nDofCon:(((size(KpropLineal,1)+1-nDofCon)+size(KUnOrd,1)-1)-nDofCon)+size(KpropRotados,1); % Dofs de KpropRotados mapeados en dofs TOTALES
dofsKtotalKSopOrd = ((((size(KpropLineal,1)+1-nDofCon)+size(KUnOrd,1)-1)-nDofCon)+size(KpropRotados,1)+1)-nDofCon:((((size(KpropLineal,1)+1-nDofCon)+size(KUnOrd,1)-1)-nDofCon)+size(KpropRotados,1)+1-nDofCon)+size(KSopOrd,1)-1; % Dofs de KSopOrd mapeados en dofs TOTALES

if  size(KpropLineal,1)-size(dofsKtotalKpropLineal,2) == 0 && size(dofsKtotalKUnOrd,2)-size(KUnOrd,1) == 0 && size(dofsKtotalKpropRotados,2)-size(KpropRotados,1) == 0 && size(dofsKtotalKSopOrd,2)-size(KSopOrd,1) == 0 && max(dofsKtotalKSopOrd)- nDofKtotal == 0
    disp('Mapeo de dofs globales correcto');
    disp('------------------------------');
else
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAPEO DE DOFS GLOBALES INCORRECTO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
end

% Sub rutina para pre - seleccionar los dof que van con BC en la malla "Alivianamientos lineales":
% ROLLER ZZ:
dofsPrim = 1:1:3*size(nodosBii,1); % De estos hay que restringir solo los (1 2 -) (4 5 -) (7 8 -) ...

% Sub rutina para pre - seleccionar los dof que van con BC en la malla "Union":
% ROLLER XX:
DofsKtotalKUnOrd = [ordenUn; dofsKtotalKUnOrd]'; % Matriz que vincula dofs propios de KUnOrd, con Ktotal
buscadorUn = ismember(ordenUn,InteriorDofsBC); % Indices de dofs propios de KUnOrd con condiciones de borde RollerXXDofsInf en ordenUn.
dofsInf = dofsKtotalKUnOrd(buscadorUn); % De estos hay que restringir solo los (- 2 3) (- 5 6) (- 8 9) ...

% Sub rutina para pre - seleccionar los dof que van con BC en la malla "Soporte":
% ROLLER ZZ:
DofsKtotalKSopOrd = [ordenSop; dofsKtotalKSopOrd ]'; % Matriz que vincula dofs propios de KUnOrd, con Ktotal
buscadorSop = ismember(ordenSop,boundaryCondSopDofs)'; % Indices de dofs propios de KUnOrd con condiciones de borde RollerXXDofsInf en ordenUn.
dofsSup = DofsKtotalKSopOrd(buscadorSop,2); % De estos hay que restringir solo los (1 2 -) (4 5 -) (7 8 -) ...

% 8b - Aplico condiciones de borde:
bcDofsInf = dofsPrim ; % Indice de Dofs globales de la punta INFERIOR
bcDofsInfPunta = dofsInf'; % Indice de Dofs globales de la punta INFERIOR derecha de la esctructura
bcDofsSup = dofsSup;% Indice de Dofs globales de la punta SUPEROPR izquierda de la esctructura

BCdofsInf = reshape(bcDofsInf,3,[])'; % Dofs reshaped en u v w columnas
BCdofsInfPunta = reshape(dofsInf,3,[])'; % Dofs reshaped en u v w columnas
BCdofsSup = reshape(bcDofsSup ,3,[])';% Dofs reshaped en u v w columnas

% Inicializo condiciones de borde:
bc = false(nDofTotal,1); 
bc([BCdofsInf(:,1) BCdofsInf(:,2)],1) = true; % RollerZZ
bc(BCdofsInfPunta(:,2),1) = true; % RollerXX
bc(BCdofsInfPunta(:,3),1) = true; % RollerXX
bc([BCdofsSup(:,1) BCdofsSup(:,2)],1) = true; % RollerZZ

%% 9 - Solver:

% 9a - Mapeo de cargas que enfrentan a Ktotal: Las de los rotados.
RsTotal = zeros(size(Ktotal,1),1);

RsTotal(dofsKtotalKpropRotados,1) = RsTotal(dofsKtotalKpropRotados,1) + Rtotal; 

% Cargas de prueba: Puntuales:
RsTotal([30 36],1) = Rpunt; % Cargas en la punta mas inferior [N]

isFree = not(bc);
isFixed = bc;
DTotal = zeros(size(Ktotal,1),1); % Son los desplazamientos de la malla COMPLETA: SE + UN + SE + AP
tic
DParcial = Ktotal(isFree,isFree)\RsTotal(isFree);

DTotal(isFree,1) = DParcial;
Dmin = min(DTotal);

disp('Dofs totales');
disp(size(Ktotal,1))

%% 10 - Post-Pro: 
   
%% 10a - Reacciones:
Rv =  Ktotal(isFixed,isFree)*DTotal(isFree);
Reacciones = zeros(nDofTotal,1);
Reacciones(isFixed) = Rv;
sumaRs = sum(RsTotal);
sumaReacciones = sum(Reacciones); 

% Vector COMPLETO DE FUERZAS:
RsCompleto = RsTotal + Reacciones; % Completo con los datos de Fuerzas que faltan

% Chequeo de reacciones:
if sumaRs+sumaReacciones < 1e-2
    disp('Reacciones bien calculadas');
    disp('--------------------------');
else
    disp('%%%%%%%%%%%%% REACCIONES NO SUMAN CERO %%%%%%%%%%%%%');
end

%% 10b.1 - Desplazamientos y tensiones sub Malla "Alivianamientos Lineales":

nDofBoundTotInf = size(KpropLineal,1); % Numero de dof boundary TOTALES de mi malla de super elementos: Viga inferior.
nDofInterior = size(interiorDofs,1); % Numero de dof interiores de UNA malla semilla.

nDofTotMeshProp = nDofBoundTotInf + NsuperElemLineal*nDofInterior; % Cantidad de dof totales de la malla de: Viga inferior.
dofsCompleteMeshLineal = zeros(nDofTotMeshProp,1); % Cantidad de desplazamientos (dof) totales (u v w)' en columna

% Mapeo de dofs con lo que mapie Ktotal  en 8a: Desplazamientos y cargas.
DboundTotalVigaInf = DTotal(dofsKtotalKpropLineal); % Desplazamientos (dof) totales (u v w)' en columna que pertenecen a esta secciòn de la malla
RboundSEInf = RsCompleto(dofsKtotalKpropLineal);
RintSE = zeros(nDofTot,1); % Vector de carga de nodos internos para cada malla semialla "Ailvianamientos"
aux1 = 1;
aux2 = nDofTot;  

% Para plotear deformada:
despX = 200*ones(size(nodos,1),1);
propX = zeros(size(nodos));
propX(:,1) = propX(:,1) + despX;

% Recuperación de tensiones en los puntos de gauss o en los nodos:
    
%    a = 1/sqrt(3); % Gauss Points
   a = 1; % Nodos
   %        ksi eta zeta
   upg = a*[-1  -1   1
            -1   1   1
             1  -1   1
             1   1   1
            -1  -1  -1
            -1   1  -1
             1  -1  -1
             1   1  -1];
      
   npg = size(upg,1);
   
% Inicializo tensione a plotear:   
avgStress = zeros(nNod,6,NsuperElemLineal);  
avgVonMises = zeros(nNod,6,NsuperElemLineal);
tic
% Para cada malla semilla:
for iMesh = 1:NsuperElemLineal
   
   % Selector de nodos Boundary:
   dofSEaB = nDofBoundaryFace*(iMesh-1) + 1;
   dofSEbB = nDofBoundaryFace*(iMesh+1); 
   dofSEB = dofSEaB:dofSEbB; % Agarro de a dos fronteras
   % Descondensación para cada malla semilla: Tomo los dofs y Rs.
   DofsBoundary = DboundTotalVigaInf(dofSEB); % Deplazamientos boundary TOMADOS de a dos fronteras: Por cada SE.
  
   DofsBoundaryI = DofsBoundary(1:nDofBoundaryFace); % Desplazamientos de los nodosBii de una malla semilla
   DofsBoundaryD = DofsBoundary(nDofBoundaryFace+1:end); % Desplazamientos de los nodosBdd de una malla semilla

   % Descondensaciòn:
   % Desplazamientos interior de cada celda malla semilla "Alivianamientos"  
   DofsInterior = (inv(K(interiorDofs,interiorDofs)))*(RintSE(interiorDofs)-K(interiorDofs,boundaryDofs)*DofsBoundary); 
   

   % Re-armado de Dof completos en columna: Para cada malla semilla.
   DofSE = [DofsBoundaryI' DofsInterior' DofsBoundaryD']';
   
   % Hago reshape de los desplazameintos para cada malla semilla para que
   % tenga la misma forma de la matriz nodos: 
   
   DboundaryIResh =  reshape(DofsBoundaryI,3,[])'; % Desp de nodosBii con su orden
   DInteriorResh  =  reshape(DofsInterior,3,[])'; % Desp de nodosI con su orden
   DboundaryDResh =  reshape(DofsBoundaryD,3,[])'; % Desp de nodosBdd con su orden
   
   % Configuraciòn de deformada para cada malla semilla "Alivianamientos" :
   % Lineales:
   
   nodePosition = zeros(size(nodos));
   
   nodePosition(nodosBii,:) = nodos(nodosBii,:) + scale*DboundaryIResh; 
   nodePosition(nodosI,:)   = nodos(nodosI,:)   + scale*DInteriorResh;
   nodePosition(nodosBdd,:) = nodos(nodosBdd,:) + scale*DboundaryDResh;
   
    % Tensiones:
   
   d = zeros(size(nodos));
   
   d(nodosBii,:) = d(nodosBii,:) + DboundaryIResh;
   d(nodosI,:)   = d(nodosI,:)   + DInteriorResh;
   d(nodosBdd,:) = d(nodosBdd,:) + DboundaryDResh;
   
   D = reshape(d',[],1);
        
   Stress = zeros(nNodEle,6,nel);
   VonMises = zeros(nNodEle,1,nel);

 for iele = 1:nel
    nodosEle = nodos(elementos(iele,:),:);
    for ipg = 1:npg
        
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        
        % Derivadas de las funciones de forma respecto de ksi, eta y zeta
         dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac  = dN*nodosEle;
        
        % Derivadas de las funciones de forma respecto de x,y,z.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        
        B(1,1:3:nDofNod*nNodEle-2) = dNxy(1,:);
        B(2,2:3:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,3:3:nDofNod*nNodEle) = dNxy(3,:);
        B(4,1:3:nDofNod*nNodEle-2) = dNxy(2,:);
        B(4,2:3:nDofNod*nNodEle-1) = dNxy(1,:);
        B(5,2:3:nDofNod*nNodEle-1) = dNxy(3,:);
        B(5,3:3:nDofNod*nNodEle) = dNxy(2,:);
        B(6,1:3:nDofNod*nNodEle-2) = dNxy(3,:);
        B(6,3:3:nDofNod*nNodEle) = dNxy(1,:);
        
        eleDoscale = nodeDoscale(elementos(iele,:),:);
        eleDoscale = reshape(eleDoscale',[],1);
        
        Stress(ipg,:,iele) = C*B*D(eleDoscale); 
        
        % Von Mises:
        SumStress = (Stress(ipg,1,iele)*Stress(ipg,2,iele) + Stress(ipg,2,iele)*Stress(ipg,3,iele) + Stress(ipg,3,iele)*Stress(ipg,1,iele));
        SumShear = (Stress(ipg,4,iele)^2 + Stress(ipg,5,iele)^2 + Stress(ipg,6,iele)^2);
        
        VonMises(ipg,:,iele) = sqrt(Stress(ipg,1,iele)^2 + Stress(ipg,2,iele)^2 + Stress(ipg,3,iele)^2 - SumStress + 3*SumShear);
    end
end
   
  % Promediado de tensiones:


for inode = 1:nNod
    [I,J] = find(elementos == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStress(inode,:,iMesh) = avgStress(inode,:,iMesh) + squeeze(Stress(J(ishare),S,I(ishare)));
        avgVonMises(inode,:,iMesh) = avgVonMises(inode,:,iMesh) + squeeze(VonMises(J(ishare),:,I(ishare)));
    end
    avgStress(inode,:,iMesh) = avgStress(inode,:,iMesh) / nShare;
    avgVonMises(inode,:,iMesh) = avgVonMises(inode,:,iMesh) / nShare;
end
 
   % Ploteo de deformadas:
   figure(4)
   Meshplot(elementos,nodos,'k',0)
   title('Malla Completa') 
   
   figure(5)
   Meshplot(elementos,nodePosition,defocolor,0)
   Meshplot(elementos,nodos,inferiorcolor,0)
   title('Deformada Completa') 
   
   figure(6)
   Meshplot(elementos,nodePosition,defocolor,0)
   Meshplot(elementos,nodos,meshcolor,0)
   title('Deformada: Sección inferior') 
   
   % Propago malla 200 mm hacia la derecha:
   nodos = nodos - propX;
   
   % Cargo y guardo los desplazamientos totales de la Viga Inferior:
   nDofSE = size(DofSE,1); % Cantidad de dofs por superelemento
   dofsCompleteMeshLineal(aux1:aux2) = dofsCompleteMeshLineal(aux1:aux2) + DofSE;
   aux1 = aux2 - nDofBoundaryFace+1;
   aux2 = aux1 + nDofTot-1;
   
end
disp('Tiempo para descondensado de viga inferior y tensiones');
toc
% Promediado de tensiones en las caras de cada SE:

% Promedio en las caras:
for iMesh = 1:NsuperElemLineal-1
   
    aux1 = (avgStress(nodosBdd,:,iMesh) + avgStress(nodosBii,:,iMesh+1))/2;
    avgStress(nodosBdd,:,iMesh) = aux1;
    avgStress(nodosBii,:,iMesh+1) = aux1;
    
    aux2 = (avgVonMises(nodosBdd,:,iMesh) + avgVonMises(nodosBii,:,iMesh+1))/2;
    avgVonMises(nodosBdd,:,iMesh) = aux2;
    avgVonMises(nodosBii,:,iMesh+1) = aux2;
 
end

% Ploteo de tenciones promediadas en las caras:

% Para plotear malla propagada:
despX = 200*ones(size(nodosPlot,1),1);
propX = zeros(size(nodosPlot));
propX(:,1) = propX(:,1) + despX;

for iMesh = 1:NsuperElemLineal
  
   tensionAVG = reshape(avgStress(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlot,elementos,tensionAVG,7);
   colormap('jet');
   axis equal
   SetColorbar 
   title('Tension XX') 

   tensionAVG = reshape(avgStress(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlot,elementos,tensionAVG,7);
   colormap('jet');
   axis equal
   SetColorbar 
   title('Tension XX') 
   
   
   tensionAVGvm = reshape(avgVonMises(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlot,elementos,tensionAVGvm,8);
   colormap('jet');
   axis equal
   SetColorbar 
   title('Tension de Von Mises en la viga inferior') 
   
   tensionAVGvm = reshape(avgVonMises(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlot,elementos,tensionAVGvm,9);
   colormap('jet');
   axis equal
   SetColorbar 
   title('Tension de Von Mises') 
   
   % Propago malla 200 mm hacia la derecha:
   nodosPlot = nodosPlot - propX;
end

% Recupero malla semilla para plotear lo que sigue:
nodos = nodos + propX; % Para recuperar las utlimas coordenadas ploteadas
nodosPlot = nodosPlot + propX; % Para recuperar las utlimas coordenadas ploteadas
clear iMesh
dofsCompleteMeshReshapedLineal = reshape(dofsCompleteMeshLineal,3,[])';

%% 10b.2 - Desplazamientos y tensiones sub Malla "Unión":
%---------------------------------------------------------------------------------------------
% NOTA: Información util:
% InferiorUnDofs = reshape(dofUn(nodosBUnInf,:)',[],1)'; % dof de los nodos inferiores ORDENADOS
% SuperiorUnDofs = reshape(dofUn(nodosBUnSup,:)',[],1)'; % dof boundary izq ORDENADOS
% InteriorDofs = reshape(dofUn(nodosUnInt,:)',[],1)'; % dof de los nodos inferiores ORDENADOS
% InteriorDofsBC = reshape(dofUn(nodosBUnbc,:)',[],1)'; % dof de los nodos inferiores ORDENADOS
% ordenUn = [InferiorUnDofs InteriorDofs InteriorDofsBC SuperiorUnDofs];
%----------------------------------------------------------------------------------------------

% Preselecciono los dofs del vector de desplazamientos totales que voy a usar:

DespUnion = DTotal(dofsKtotalKUnOrd); % Desplazamientos de la Union en idioma Ktotal.
DespUnResh = reshape(DespUnion,3,[])';
% Re-mapeo de los dofs en idioma Un: En columna y rspetando el ordenUn.
InfConDofs = 1:size(InferiorUnDofs,2);
IntDofs = size(InferiorUnDofs,2)+1:size(InferiorUnDofs,2)+size(InteriorDofs,2);
IntDofsBC = (size(InferiorUnDofs,2)+size(InteriorDofs,2))+1:(size(InferiorUnDofs,2)+size(InteriorDofs,2))+size(InteriorDofsBC,2);
SupConDofs = ((size(InferiorUnDofs,2)+size(InteriorDofs,2))+size(InteriorDofsBC,2))+1:((size(InferiorUnDofs,2)+size(InteriorDofs,2))+size(InteriorDofsBC,2))+size(SuperiorUnDofs,2);

DofsBoundInfUnCOL = DespUnion(InfConDofs,1); 
DofsInteriorDofsUnCOL = DespUnion(IntDofs,1);
DofsBoundBCUnCOL = DespUnion(IntDofsBC,1);
DofsBoundSupUnCOL = DespUnion(SupConDofs,1);

% Reshape de los desplazamientos (u v w): Para hablar en idioma Nodos.
DofsBoundInfUn = reshape(DofsBoundInfUnCOL,3,[])';
DofsInteriorDofsUn = reshape(DofsInteriorDofsUnCOL,3,[])';
DofsBoundBCUn = reshape(DofsBoundBCUnCOL,3,[])';
DofsBoundSupUn = reshape(DofsBoundSupUnCOL,3,[])';

% Configuracion de deformada:
nodePositionUn = zeros(size(nodosUn));

% Para plotear deformada:
nodALIV = nodos(nodosBdd(1,:),1); % Ultimas coordenadas X de nodo de la esquina inf DERECHA de ALIVIANAMIENTOS Lineal 
nodUn = nodosUn(nodosBUnInf(1,:),1); % Idem perod de Union de la parte de abajo.
difUn = nodALIV-nodUn;
despXUn = difUn*ones(size(nodosUn,1),1);
propXUn = zeros(size(nodosUn));
propXUn(:,1) = propXUn(:,1) + despXUn;

%-----------------------------------------------------------------------
% NOTA: La malla "ALIVIANAMIENTOS" Lineal se dedujo con el (0,0,0) de la
% misma en el centro de la malla, para facilitar la sub rutina de
% condensacon.
%-----------------------------------------------------------------------

nodosUn = nodosUn + propXUn; % Acerco la malla Union para empalmar geometricamente con ALIVIANAMIENTOS LINEAL

nodePositionUn(nodosBUnInf,:) = nodosUn(nodosBUnInf,:)  + scale*DofsBoundInfUn; %DboundaryDResh
nodePositionUn(nodosUnInt,:)  = nodosUn(nodosUnInt,:)   + scale*DofsInteriorDofsUn;
nodePositionUn(nodosBUnbc,:)  = nodosUn(nodosBUnbc,:)   + scale*DofsBoundBCUn;
nodePositionUn(nodosBUnSup,:) = nodosUn(nodosBUnSup,:)  + scale*DofsBoundSupUn;

% Ploteo de deformadas:
figure(4)
Meshplot(elementosUn,nodosUn,'k',0)

figure(5)
Meshplot(elementosUn,nodePositionUn,defocolor,0)
Meshplot(elementosUn,nodosUn,unioncolor,0)

figure(10)
Meshplot(elementosUn,nodePositionUn,defocolor,0)
Meshplot(elementosUn,nodosUn,meshcolor,0)
title('Deformada: Sección unión') 

% Tensiones:
   
   dUn = zeros(size(nodosUn));
   
   dUn(nodosBUnInf,:)  = dUn(nodosBUnInf,:) + DofsBoundInfUn;
   dUn(nodosUnInt,:)   = dUn(nodosUnInt,:)  + DofsInteriorDofsUn;
   dUn(nodosBUnbc,:)   = dUn(nodosBUnbc,:)  + DofsBoundBCUn;
   dUn(nodosBUnSup,:)  = dUn(nodosBUnSup,:) + DofsBoundSupUn;
   
   DUn = reshape(dUn',[],1);
        
   StressUn = zeros(nNodEleUn,6,nelUn);
   VonMisesUn = zeros(nNodEleUn,1,nelUn);
   
for iele = 1:nelUn
    nodosEle = nodosUn(elementosUn(iele,:),:);
    for inode = 1:nNodEleUn
        
        % Punto de Gauss
        ksi = upg(inode,1);
        eta = upg(inode,2);
        zeta = upg(inode,3);
        
        % Derivadas de las funciones de forma respecto de ksi, eta y zeta
         dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac  = dN*nodosEle;
        
        % Derivadas de las funciones de forma respecto de x,y,z.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        
        B = zeros(size(C,2),nDofNodUn*nNodEleUn);
        
        B(1,1:3:nDofNodUn*nNodEleUn-2) = dNxy(1,:);
        B(2,2:3:nDofNodUn*nNodEleUn-1) = dNxy(2,:);
        B(3,3:3:nDofNodUn*nNodEleUn) = dNxy(3,:);
        B(4,1:3:nDofNodUn*nNodEleUn-2) = dNxy(2,:);
        B(4,2:3:nDofNodUn*nNodEleUn-1) = dNxy(1,:);
        B(5,2:3:nDofNodUn*nNodEleUn-1) = dNxy(3,:);
        B(5,3:3:nDofNodUn*nNodEleUn) = dNxy(2,:);
        B(6,1:3:nDofNodUn*nNodEleUn-2) = dNxy(3,:);
        B(6,3:3:nDofNodUn*nNodEleUn) = dNxy(1,:);
        
        eleDoscale = nodeDoscaleUn(elementosUn(iele,:),:);
        eleDoscale = reshape(eleDoscale',[],1);
        StressUn(inode,:,iele) = C*B*DUn(eleDoscale);
        
        % Von Mises:
        SumStressUn = (StressUn(ipg,1,iele)*StressUn(ipg,2,iele) + StressUn(ipg,2,iele)*StressUn(ipg,3,iele) + StressUn(ipg,3,iele)*StressUn(ipg,1,iele));
        SumShearUn = (StressUn(ipg,4,iele)^2 + StressUn(ipg,5,iele)^2 + StressUn(ipg,6,iele)^2);
        
        VonMisesUn(ipg,:,iele) = sqrt(StressUn(ipg,1,iele)^2 + StressUn(ipg,2,iele)^2 + StressUn(ipg,3,iele)^2 - SumStressUn + 3*SumShearUn);
    end
end

% Promediado de tensiones:
avgStressUn = zeros(nNodUn,nelUn);
avgVonMisesUn = zeros(nNodUn,nelUn);

for inode = 1:nNodUn
    [I,J] = find(elementosUn == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStressUn(inode,:) = avgStressUn(inode,:) + squeeze(StressUn(J(ishare),S,I(ishare)));
        avgVonMisesUn(inode,:) = avgVonMisesUn(inode,:) + squeeze(VonMisesUn(J(ishare),:,I(ishare)));
    end
    
    avgStressUn(inode,:) = avgStressUn(inode,:) / nShare;
    avgVonMisesUn(inode,:) = avgVonMisesUn(inode,:) /nShare;
end

% Ploteo de tension:

tensionAVGUn = reshape(avgStressUn(:,:),1,[])';
PlotFieldonMesh(nodosUn,elementosUn,tensionAVGUn,7);
colormap('jet');
axis equal
SetColorbar 

tensionAVGUn = reshape(avgStressUn(:,:),1,[])';
PlotFieldonMesh(nodosUn,elementosUn,tensionAVGUn,11);
colormap('jet');
title('Tension XX en la Union');
axis equal
SetColorbar 

tensionVmAVGUn = reshape(avgVonMisesUn(:,:),1,[])';
PlotFieldonMesh(nodosUn,elementosUn,tensionVmAVGUn ,9);
colormap('jet');
axis equal
SetColorbar 

tensionVmAVGUn = reshape(avgVonMisesUn(:,:),1,[])';
PlotFieldonMesh(nodosUn,elementosUn,tensionVmAVGUn ,12);
colormap('jet');
title('Tension de Von Mises en la Union');
axis equal
SetColorbar 

%% 10b.3 - Desplazamientos y tensiones sub Malla "Alivianamientos rotados":

% Mapeo de dofs con lo que mapie Ktotal  en 8a: Desplazamientos y cargas.
DboundTotalVigaSup = DTotal(dofsKtotalKpropRotados); % Desplazamientos (dof) totales (u v w)' en columna que pertenecen a esta secciòn de la malla
RboundSESup = RsTotal(dofsKtotalKpropRotados);
RintSESup = R; % Vector de carga de nodos internos para cada malla semialla "Ailvianamientos"
aux1 = 1;
aux2 = nDofTot;  

% Rotacion de nodos de "sub malla alivianamientos" para obtener las coordenadas de primer malla semilla rotada:

nodosRota = nodos*Lamda'; % Rotacion de matriz nodos.
difer = nodosRota(nodosBdd,:)- nodosUn(nodosBUnSup,:); % Calculo de diferencia para conectar nodos de geometría
% Traslación de matriz nodos:
diferX = difer(1,1)*ones(size(nodos,1),1);
diferZ = difer(1,3)*ones(size(nodos,1),1);

nodosRota(:,1) = nodosRota(:,1) - diferX;
nodosRota(:,3) = nodosRota(:,3) - diferZ;

% Para plotear deformada y tensiones:
% Malla para propagar:
nodosRot = nodosRota; % Matriz de coord de nodos rotados para la primer malla "alivianamientos rotados"
elementosRot = elementosOrd;
nodosPlotRot = nodosRota;

% Inicio de componentes para trasladar la malla semilla "rotada": 12 veces.
traslacionX = abs((200*cosd(theta)))*ones(size(nodos,1),1);
traslacionZ = abs((200*sind(theta)))*ones(size(nodos,1),1);

% Inicializo tensione a plotear:   
avgStressRot = zeros(nNod,6,NsuperElemRotados);  
avgVonMisesRot = zeros(nNod,6,NsuperElemRotados);
tic
for iMesh = 1:NsuperElemRotados
    
   % Selector de nodos Boundary:
   dofSEaB = nDofBoundaryFace*(iMesh-1) + 1;
   dofSEbB = nDofBoundaryFace*(iMesh+1); 
   dofSEB = dofSEaB:dofSEbB; % Agarro de a dos fronteras
   % Descondensación para cada malla semilla: Tomo los dofs y Rs.
   DofsBoundary = DboundTotalVigaSup(dofSEB); % Deplazamientos boundary TOMADOS de a dos fronteras: Por cada SE.
   
   DofsBoundaryD = DofsBoundary(1:nDofBoundaryFace); % Desplazamientos de los nodosBii de una malla semilla rotada
   DofsBoundaryI = DofsBoundary(nDofBoundaryFace+1:end); % Desplazamientos de los nodosBdd de una malla semilla rotada
   
   % Descondensaciòn:
   % Desplazamientos interior de cada celda malla semilla "Alivianamientos"  
   DofsInterior = (inv(Krot(interiorDofs,interiorDofs)))*(RintSESup(interiorDofs)-Krot(interiorDofs,boundaryDofs)*DofsBoundary); 
   
   % Re-armado de Dof completos en columna: Para cada malla semilla.
   DofSE = [DofsBoundaryI' DofsInterior' DofsBoundaryD']';
   
   % Hago reshape de los desplazameintos para cada malla semilla para que
   % tenga la misma forma de la matriz nodos: 
   
   DboundaryIResh =  reshape(DofsBoundaryI,3,[])'; % Desp de nodosBii con su orden
   DInteriorResh  =  reshape(DofsInterior,3,[])'; % Desp de nodosI con su orden
   DboundaryDResh =  reshape(DofsBoundaryD,3,[])'; % Desp de nodosBdd con su orden
   
   % Configuraciòn de deformada para cada malla semilla "Alivianamientos" :
   % Rotados:
   
   nodePositionRot = zeros(size(nodosRot));    
   
   nodePositionRot(nodosBii,:) = nodosRot(nodosBii,:) + scale*DboundaryIResh; 
   nodePositionRot(nodosI,:)   = nodosRot(nodosI,:)   + scale*DInteriorResh;
   nodePositionRot(nodosBdd,:) = nodosRot(nodosBdd,:) + scale*DboundaryDResh; 
   
   % Tensiones:
   
   dRot = zeros(size(nodosRota));
   
   dRot(nodosBii,:) = dRot(nodosBii,:) + DboundaryIResh;
   dRot(nodosI,:)   = dRot(nodosI,:)   + DInteriorResh;
   dRot(nodosBdd,:) = dRot(nodosBdd,:) + DboundaryDResh;
   
   DRot = reshape(dRot',[],1);
        
   StressRot = zeros(nNodEle,6,nel);
   VonMisesRot = zeros(nNodEle,6,nel);
   
for iele = 1:nel
    nodosEle = nodosRot(elementosRot(iele,:),:);
    for ipg = 1:npg
        
        % Punto de Gauss
        ksi = upg(ipg,1);
        eta = upg(ipg,2);
        zeta = upg(ipg,3);
        
        % Derivadas de las funciones de forma respecto de ksi, eta y zeta
         dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac  = dN*nodosEle;
        
        % Derivadas de las funciones de forma respecto de x,y,z.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        
        B = zeros(size(C,2),nDofNod*nNodEle);
        
        B(1,1:3:nDofNod*nNodEle-2) = dNxy(1,:);
        B(2,2:3:nDofNod*nNodEle-1) = dNxy(2,:);
        B(3,3:3:nDofNod*nNodEle) = dNxy(3,:);
        B(4,1:3:nDofNod*nNodEle-2) = dNxy(2,:);
        B(4,2:3:nDofNod*nNodEle-1) = dNxy(1,:);
        B(5,2:3:nDofNod*nNodEle-1) = dNxy(3,:);
        B(5,3:3:nDofNod*nNodEle) = dNxy(2,:);
        B(6,1:3:nDofNod*nNodEle-2) = dNxy(3,:);
        B(6,3:3:nDofNod*nNodEle) = dNxy(1,:);
        
        eleDoscale = nodeDoscale(elementosRot(iele,:),:);
        eleDoscale = reshape(eleDoscale',[],1);
        
         StressRot (ipg,:,iele) = C*B*DRot(eleDoscale); 
        
        % Von Mises:
        SumStressRot = (StressRot(ipg,1,iele)* StressRot(ipg,2,iele) +  StressRot(ipg,2,iele)*StressRot(ipg,3,iele) + StressRot(ipg,3,iele)*StressRot(ipg,1,iele));
        SumShearRot =  (StressRot(ipg,4,iele)^2 + StressRot(ipg,5,iele)^2 + StressRot(ipg,6,iele)^2);
        
        VonMisesRot(ipg,:,iele) = sqrt(StressRot(ipg,1,iele)^2 + StressRot(ipg,2,iele)^2 + StressRot(ipg,3,iele)^2 - SumStressRot + 3*SumShearRot);
    end
end
   
% Promediado de tensiones:


for inode = 1:nNod
    [I,J] = find(elementosRot == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStressRot(inode,:,iMesh) = avgStressRot(inode,:,iMesh) + squeeze(StressRot(J(ishare),S,I(ishare)));
        avgVonMisesRot(inode,:,iMesh) = avgVonMisesRot(inode,:,iMesh) + squeeze(VonMisesRot(J(ishare),S,I(ishare)));
    end
    avgStressRot(inode,:,iMesh) = avgStressRot(inode,:,iMesh) / nShare;
    avgVonMisesRot(inode,:,iMesh) = avgVonMisesRot(inode,:,iMesh) / nShare;
end
   
   % Ploteo de deformada:
   figure(4)
   Meshplot(elementosRot,nodosRot,'k',0)

   figure(5)
   Meshplot(elementosRot,nodePositionRot,defocolor,0)
   Meshplot(elementosRot,nodosRot,superiorcolor,0)
   
   figure(13)
   Meshplot(elementosRot,nodePositionRot,defocolor,0)
   Meshplot(elementosRot,nodosRot,meshcolor,0)
   title('Deformada: Sección superior') 
   
   % Desplazo malla hacia arriba a la izquierda:

   nodosRot(:,1) = nodosRot(:,1) + traslacionX;
   nodosRot(:,3) = nodosRot(:,3) + traslacionZ;


end
disp('Tiempo de descondensado de malla superior con tensiones');
toc
% Promediado de tensiones en las caras: 

% Teniendo en cuenta la cara superior de la unión:
avgStressBUnSup = avgStressUn(nodosBUnSup,35:40);
avgVonMisesBUnSup =  avgVonMisesUn(nodosBUnSup,35:40);

avgStressRot(nodosBdd,:,1)   = (avgStressRot(nodosBdd,:,1) + avgStressBUnSup)/2;
avgVonMisesRot(nodosBdd,:,1) = (avgVonMisesRot(nodosBdd,:,1) + avgVonMisesBUnSup)/2;

for iMesh = 1:NsuperElemRotados-1
   
    aux1 = (avgStressRot(nodosBdd,:,iMesh) + avgStressRot(nodosBii,:,iMesh+1))/2;
    avgStressRot(nodosBdd,:,iMesh) = aux1;
    avgStressRot(nodosBii,:,iMesh+1) = aux1;
    
    aux2 = (avgVonMisesRot(nodosBdd,:,iMesh) + avgVonMisesRot(nodosBii,:,iMesh+1))/2;
    avgVonMisesRot(nodosBdd,:,iMesh) = aux2;
    avgVonMisesRot(nodosBii,:,iMesh+1) = aux2;
 
end

% Ploteo de tenciones promediadas en las caras de los SE.

% Inicio de componentes para trasladar la malla semilla "rotada": 12 veces.
traslacionX = abs((200*cosd(theta)))*ones(size(nodos,1),1);
traslacionZ = abs((200*sind(theta)))*ones(size(nodos,1),1);


for iMesh = 1:NsuperElemRotados
  
    
   tensionAVGRot = reshape(avgStressRot(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlotRot,elementosRot,tensionAVGRot,7);
   colormap('jet')
   axis equal
   SetColorbar 
   
   tensionAVGvmRot = reshape(avgVonMisesRot(:,:,iMesh),1,[])';
   PlotFieldonMesh(nodosPlotRot,elementosRot,tensionAVGvmRot,9);
   colormap('jet')
   axis equal
   SetColorbar 
   

      % Desplazo malla hacia arriba a la izquierda:
  
   nodosPlotRot(:,1) = nodosPlotRot(:,1) + traslacionX;
   nodosPlotRot(:,3) = nodosPlotRot(:,3) + traslacionZ;

end

% Re acomodo el ultimo super elemento para plotear "Soporte":
nodosRot(:,1) = nodosRot(:,1) - traslacionX;
nodosRot(:,3) = nodosRot(:,3) - traslacionZ;

%% 10b.4 - Desplazamientos y tensiones sub Malla "Soporte":
% Selecciono los desplazamientos:
DespSop = DTotal(dofsKtotalKSopOrd,1);
diferS = nodosRot(nodosBii(1,:),:) - nodosSop(nodosConectoresSop(1,:),:); % Calculo de diferencia para conectar nodos de geometría
% Traslación de matriz nodos:
diferXS = diferS(1,1)*ones(size(nodosSop,1),1);
diferZS = diferS(1,3)*ones(size(nodosSop,1),1);

nodosSop(:,1) = nodosSop(:,1) + diferXS;
nodosSop(:,3) = nodosSop(:,3) + diferZS;

nodosSopMesh = nodosSop; % Matriz de coord de nodos de la geometría global

%---------------------------------------------------------------------------------------------------------------------------------------
%NOTA: Info util
% conectorSopDofsOrd = reshape(dofSop(nodosConectoresSop,:)',[],1)'; % dof que conectan con la viga superior Ordenados criterio IZQUIERDA
% boundaryCondSopDofs = reshape(dofSop(nodosConBCSop,:)',[],1)'; % dof que toman BC
% interiorSopDofs = reshape(dofSop(nodosInternosSop,:)',[],1)'; % dof de los nodos interiores
% ordenSop = [conectorSopDofsOrd interiorSopDofs boundaryCondSopDofs];
%-----------------------------------------------------------------------------------------------------------------------------------------

% Re - Mapeo de los desplazamientos:
DofsCon = 1:size(conectorSopDofsOrd,2);
IntDofs = size(conectorSopDofsOrd,2)+1:size(conectorSopDofsOrd,2)+size(interiorSopDofs,2);
BoundBCDofs = (size(conectorSopDofsOrd,2)+size(interiorSopDofs,2))+1:(size(conectorSopDofsOrd,2)+size(interiorSopDofs,2))+size(boundaryCondSopDofs,2);

DespConCOL = DespSop(DofsCon,1);
DespIntCOL = DespSop(IntDofs,1);
DespBCCOL = DespSop(BoundBCDofs,1);

DespConResh = reshape(DespConCOL,3,[])';
DespIntResh = reshape(DespIntCOL,3,[])';
DespBCResh =  reshape(DespBCCOL,3,[])';

nodePositionSop = zeros(size(nodosSop));

nodePositionSop(nodosConectoresSop,:) = nodosSopMesh(nodosConectoresSop,:)  + scale*DespConResh;
nodePositionSop(nodosInternosSop,:)   = nodosSopMesh(nodosInternosSop,:)    + scale*DespIntResh;
nodePositionSop(nodosConBCSop,:)      = nodosSopMesh(nodosConBCSop,:)       + scale*DespBCResh;

% Tensiones:

dSop = zeros(size(nodosSop));

dSop(nodosConectoresSop,:) = dSop(nodosConectoresSop,:) + DespConResh;
dSop(nodosInternosSop,:)   = dSop(nodosInternosSop,:)   + DespIntResh;
dSop(nodosConBCSop,:)      = dSop(nodosConBCSop,:)      + DespBCResh;

DSop = reshape(dSop',[],1);
        
StressSop = zeros(nNodEleSop,6,nelSop);
VonMisesSop = zeros(nNodEleSop,1,nelSop);
   
for iele = 1:nelSop
    nodosEle = nodosSop(elementosSop(iele,:),:);
    for inode = 1:nNodEleSop
        
        % Punto de Gauss
        ksi = upg(inode,1);
        eta = upg(inode,2);
        zeta = upg(inode,3);
        
        % Derivadas de las funciones de forma respecto de ksi, eta y zeta
         dN = shapefunsder_3D([ksi eta zeta],eleType);
        
        % Derivadas de x,y,z respecto de ksi, eta y zeta
        jac  = dN*nodosEle;
        
        % Derivadas de las funciones de forma respecto de x,y,z.
        dNxy = jac\dN;          % dNxy = inv(jac)*dN
        
        % Matriz de derivadas B: Para 3D
        
        B = zeros(size(C,2),nDofNodSop*nNodEleSop);
        
        B(1,1:3:nDofNodSop*nNodEleSop-2) = dNxy(1,:);
        B(2,2:3:nDofNodSop*nNodEleSop-1) = dNxy(2,:);
        B(3,3:3:nDofNodSop*nNodEleSop) = dNxy(3,:);
        B(4,1:3:nDofNodSop*nNodEleSop-2) = dNxy(2,:);
        B(4,2:3:nDofNodSop*nNodEleSop-1) = dNxy(1,:);
        B(5,2:3:nDofNodSop*nNodEleSop-1) = dNxy(3,:);
        B(5,3:3:nDofNodSop*nNodEleSop) = dNxy(2,:);
        B(6,1:3:nDofNodSop*nNodEleSop-2) = dNxy(3,:);
        B(6,3:3:nDofNodSop*nNodEleSop) = dNxy(1,:);
        
        eleDoscale = nodeDoscaleUn(elementosSop(iele,:),:);
        eleDoscale = reshape(eleDoscale',[],1);
        StressSop(inode,:,iele) = C*B*DSop(eleDoscale);
        
        % Von Mises:
        SumStressSop = (StressSop(ipg,1,iele)*StressSop(ipg,2,iele) + StressSop(ipg,2,iele)*StressSop(ipg,3,iele) + StressSop(ipg,3,iele)*StressSop(ipg,1,iele));
        SumShearSop = (StressSop(ipg,4,iele)^2 + StressSop(ipg,5,iele)^2 + StressSop(ipg,6,iele)^2);
        
        VonMisesSop(ipg,:,iele) = sqrt(StressSop(ipg,1,iele)^2 + StressSop(ipg,2,iele)^2 + StressSop(ipg,3,iele)^2 - SumStressSop + 3*SumShearSop);
    end
end

% Promediado de tensiones:
avgStressSop = zeros(nNodSop,nelSop);
avgVonMisesSop = zeros(nNodSop,nelSop);

for inode = 1:nNodSop
    [I,J] = find(elementosSop == inode);
    nShare = length(I);
    for ishare = 1:nShare
        avgStressSop(inode,:) = avgStressSop(inode,:) + squeeze(StressSop(J(ishare),S,I(ishare)));
        avgVonMisesSop(inode,:) = avgVonMisesSop(inode,:) + squeeze(VonMisesSop(J(ishare),:,I(ishare)));
    end
    
    avgStressSop(inode,:) = avgStressSop(inode,:) / nShare;
    avgVonMisesSop(inode,:) = avgVonMisesSop(inode,:) /nShare;
end

% Ploteo de tension:

tensionAVGSop = reshape(avgStressSop(:,:),1,[])';
PlotFieldonMesh(nodosSop,elementosSop',tensionAVGSop,7);
colormap('jet');
axis equal
SetColorbar 

tensionAVGSop = reshape(avgStressSop(:,:),1,[])';
PlotFieldonMesh(nodosSop,elementosSop',tensionAVGSop,11);
colormap('jet');
title('Tension XX en el soporte');
axis equal
SetColorbar 

tensionVmAVGSop = reshape(avgVonMisesSop(:,:),1,[])';
PlotFieldonMesh(nodosSop,elementosSop',tensionVmAVGSop ,9);
colormap('jet');
axis equal
SetColorbar 

tensionVmAVGSop = reshape(avgVonMisesSop(:,:),1,[])';
PlotFieldonMesh(nodosSop,elementosSop',tensionVmAVGSop ,12);
colormap('jet');
title('Tension de Von Mises: Malla completa');
axis equal
SetColorbar 

% Ploteo de deformada:
figure(4)
Meshplot(elementosSop,nodosSop,'k',0)

figure(5)
Meshplot(elementosSop,nodePositionSop,defocolor,0)
Meshplot(elementosSop,nodosSop,sopcolor,0)

figure(14)
Meshplot(elementosSop,nodePositionSop,defocolor,0)
Meshplot(elementosSop,nodosSop,meshcolor,0)
title('Deformada: Sección soporte') 

toc