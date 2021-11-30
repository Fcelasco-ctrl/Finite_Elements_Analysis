%  PARA  PROBLEMA 3D
function [Ni,N,Nf] = shapefuns_3D(pointArray,eleType)


ngauss = size(pointArray,1);

switch eleType

    case 'Q9'
        N  = zeros(2,18,ngauss);
        Ni = zeros(1, 9,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N9 =      (1 - ksi^2)*(1 - eta^2);
            N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
            N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
            N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
            N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);

            Ni(1,:     ,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (1,1:2:17,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (2,2:2:18,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
        end

    case 'AHMAD9'
        N  = zeros(3,3*9,ngauss);
        Ni = zeros(1,  9,ngauss);
        Id = eye(3);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N9 =      (1 - ksi^2)*(1 - eta^2);
            N8 = 0.50*(1 - ksi  )*(1 - eta^2) - 0.5*N9;
            N7 = 0.50*(1 - ksi^2)*(1 + eta  ) - 0.5*N9;
            N6 = 0.50*(1 + ksi  )*(1 - eta^2) - 0.5*N9;
            N5 = 0.50*(1 - ksi^2)*(1 - eta  ) - 0.5*N9;
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8 + 0.5*N9);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7 + 0.5*N9);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6 + 0.5*N9);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8 + 0.5*N9);

            Ni(1,:,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8 N9];
            N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, ...
                N4*Id, N5*Id, N6*Id, ...
                N7*Id, N8*Id, N9*Id ];
        end

    case 'Q8'
        N  = zeros(2,16,ngauss);
        Ni = zeros(1, 8,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N8 = 0.50*(1 - ksi  )*(1 - eta^2);
            N7 = 0.50*(1 - ksi^2)*(1 + eta  );
            N6 = 0.50*(1 + ksi  )*(1 - eta^2);
            N5 = 0.50*(1 - ksi^2)*(1 - eta  );
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);

            Ni(1,:     ,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (1,1:2:15,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (2,2:2:16,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
        end

    case 'AHMAD8'
        N  = zeros(3,3*8,ngauss);
        Ni = zeros(1,  8,ngauss);
        Id = eye(3);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N8 = 0.50*(1 - ksi  )*(1 - eta^2);
            N7 = 0.50*(1 - ksi^2)*(1 + eta  );
            N6 = 0.50*(1 + ksi  )*(1 - eta^2);
            N5 = 0.50*(1 - ksi^2)*(1 - eta  );
            N4 = 0.25*(1 - ksi  )*(1 + eta  ) - 0.5*(N7 + N8);
            N3 = 0.25*(1 + ksi  )*(1 + eta  ) - 0.5*(N6 + N7);
            N2 = 0.25*(1 + ksi  )*(1 - eta  ) - 0.5*(N5 + N6);
            N1 = 0.25*(1 - ksi  )*(1 - eta  ) - 0.5*(N5 + N8);

            Ni(1,:,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
            N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, ...
                N4*Id, N5*Id, N6*Id, ...
                N7*Id, N8*Id ];
        end

    case 'Q4'
        N  = zeros(2,8,ngauss);
        Ni = zeros(1,4,ngauss);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N4 = 0.25*(1 - ksi)*(1 + eta);
            N3 = 0.25*(1 + ksi)*(1 + eta);
            N2 = 0.25*(1 + ksi)*(1 - eta);
            N1 = 0.25*(1 - ksi)*(1 - eta);

            Ni(1,:    ,igauss) = [N1 N2 N3 N4];
            N (1,1:2:7,igauss) = [N1 N2 N3 N4];
            N (2,2:2:8,igauss) = [N1 N2 N3 N4];
        end
        
    case 'H8'  % Segun Cook pag 217.
         N  = zeros(3,24,ngauss);
         Ni = zeros(1,8,ngauss);
         for igauss = 1:ngauss
             ksi  = pointArray(igauss,1);
             eta  = pointArray(igauss,2);
             zeta = pointArray(igauss,3);
             
             N8 = eta/8 + ksi/8 + zeta/8 + (eta*ksi)/8 + (eta*zeta)/8 + (ksi*zeta)/8 + (eta*ksi*zeta)/8 + 1/8;
             N7 = eta/8 + ksi/8 - zeta/8 + (eta*ksi)/8 - (eta*zeta)/8 - (ksi*zeta)/8 - (eta*ksi*zeta)/8 + 1/8;
             N6 = ksi/8 - eta/8 - zeta/8 - (eta*ksi)/8 + (eta*zeta)/8 - (ksi*zeta)/8 + (eta*ksi*zeta)/8 + 1/8;
             N5 = ksi/8 - eta/8 + zeta/8 - (eta*ksi)/8 - (eta*zeta)/8 + (ksi*zeta)/8 - (eta*ksi*zeta)/8 + 1/8;
             N4 = eta/8 - ksi/8 + zeta/8 - (eta*ksi)/8 + (eta*zeta)/8 - (ksi*zeta)/8 - (eta*ksi*zeta)/8 + 1/8;
             N3 = eta/8 - ksi/8 - zeta/8 - (eta*ksi)/8 - (eta*zeta)/8 + (ksi*zeta)/8 + (eta*ksi*zeta)/8 + 1/8;
             N2 = (eta*ksi)/8 - ksi/8 - zeta/8 - eta/8 + (eta*zeta)/8 + (ksi*zeta)/8 - (eta*ksi*zeta)/8 + 1/8;
             N1 = zeta/8 - ksi/8 - eta/8 + (eta*ksi)/8 - (eta*zeta)/8 - (ksi*zeta)/8 + (eta*ksi*zeta)/8 + 1/8;
             
             Ni(1,:    ,igauss)  = [N1 N2 N3 N4 N5 N6 N7 N8];
             Nf = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (1,1:3:22,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (2,2:3:23,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (3,3:3:24,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
         end
        
         case 'H8ADINA'  % Segun manual de ADINA.
         N  = zeros(3,24,ngauss);
         Ni = zeros(1,8,ngauss);
         for igauss = 1:ngauss
             ksi  = pointArray(igauss,1);
             eta  = pointArray(igauss,2);
             zeta = pointArray(igauss,3);
             
             N8 = (1/8)*(1+ksi)*(1+eta)*(1-zeta);
             N7 = (1/8)*(ksi-1)*(1+eta)*(zeta-1);
             N6 = (1/8)*(1-ksi)*(1+eta)*(1+zeta);
             N5 = (1/8)*(1+ksi)*(1+eta)*(1+zeta);
             N4 = (1/8)*(1+ksi)*(eta-1)*(zeta-1);
             N3 = (1/8)*(1-ksi)*(eta-1)*(zeta-1);
             N2 = (1/8)*(ksi-1)*(eta-1)*(1+zeta);
             N1 = (1/8)*(1+ksi)*(1-eta)*(1+zeta);
             
             Ni(1,:    ,igauss)  = [N1 N2 N3 N4 N5 N6 N7 N8];
             Nf = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (1,1:3:22,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (2,2:3:23,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (3,3:3:24,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             
         end
         
         case 'H8FACU'  % Segun Criterio de Elements_sorter
         N  = zeros(3,24,ngauss);
         Ni = zeros(1,8,ngauss);
         for igauss = 1:ngauss
             ksi  = pointArray(igauss,1);
             eta  = pointArray(igauss,2);
             zeta = pointArray(igauss,3);
             
             N8 = (1/8)*(1+ksi)*(1+eta)*(1-zeta);
             N7 = (1/8)*(1+ksi)*(eta-1)*(zeta-1);
             N6 = (1/8)*(1-ksi)*(eta-1)*(zeta-1);
             N5 = (1/8)*(ksi-1)*(1+eta)*(zeta-1);
             N4 = (1/8)*(1+ksi)*(1+eta)*(1+zeta);
             N3 = (1/8)*(1+ksi)*(1-eta)*(1+zeta);
             N2 = (1/8)*(ksi-1)*(eta-1)*(1+zeta); 
             N1 = (1/8)*(1-ksi)*(1+eta)*(1+zeta);
             
             Ni(1,:    ,igauss)  = [N1 N2 N3 N4 N5 N6 N7 N8];
             Nf = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (1,1:3:22,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (2,2:3:23,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             N (3,3:3:24,igauss) = [N1 N2 N3 N4 N5 N6 N7 N8];
             
         end
         
    case 'AHMAD4'
        N  = zeros(3,3*4,ngauss);
        Ni = zeros(1,  4,ngauss);
        Id = eye(3);
        for igauss = 1:ngauss
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            N4 = 0.25*(1 - ksi)*(1 + eta);
            N3 = 0.25*(1 + ksi)*(1 + eta);
            N2 = 0.25*(1 + ksi)*(1 - eta);
            N1 = 0.25*(1 - ksi)*(1 - eta);

            Ni(1,:,igauss) = [N1 N2 N3 N4];
            N (:,:,igauss) = [ N1*Id, N2*Id, N3*Id, N4*Id ];
        end

end


