% PARA PROBLEMA 3D
function dN = shapefunsder_3D(pointArray,eleType)

ngauss = size(pointArray,1);

switch eleType

    case {'Q9', 'AHMAD9'}
        dN = zeros(2,9,ngauss);
        for igauss = 1:ngauss
            
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            dN(:,:,igauss) = [ % derivadas respecto de ksi
                0.25*eta*(-1+eta)*(2*ksi-1),      0.25*eta*(-1+eta)*(2*ksi+1),       0.25*eta*(1+eta)*(2*ksi+1),...
                0.25*eta*( 1+eta)*(2*ksi-1),                -ksi*eta*(-1+eta),  -1/2*(-1+eta)*(1+eta)*(2*ksi+1),...
                           -ksi*eta*(1+eta),  -1/2*(-1+eta)*(1+eta)*(2*ksi-1),           2*ksi*(-1+eta)*(1+eta)
                % derivadas respecto de eta
                   0.25*ksi*(-1+2*eta)*(ksi-1),      0.25*ksi*(-1+2*eta)*(1+ksi),       0.25*ksi*(2*eta+1)*(1+ksi),...
                    0.25*ksi*(2*eta+1)*(ksi-1),  -0.5*(ksi-1)*(1+ksi)*(-1+2*eta),                 -ksi*eta*(1+ksi),...
                -0.5*(ksi-1)*(1+ksi)*(2*eta+1),                 -ksi*eta*(ksi-1),           2*(ksi-1)*(1+ksi)*eta ];
        end

    case {'Q8', 'AHMAD8'}
        dN = zeros(2,8,ngauss);
        for igauss = 1:ngauss
            
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            dN(:,:,igauss) = [  % derivadas respecto de ksi
                -0.25*(-1+eta)*(eta+2*ksi),  -0.25*(-1+eta)*(-eta+2*ksi),    0.25*(1+eta)*(eta+2*ksi),   0.25*(1+eta)*(-eta+2*ksi),...
                              ksi*(-1+eta),        -0.5*(-1+eta)*(1+eta),                -ksi*(1+eta),        0.5*(-1+eta)*(1+eta)
                % derivadas respecto de eta
                -0.25*(-1+ksi)*(ksi+2*eta),   -0.25*(1+ksi)*(ksi-2*eta),    0.25*(1+ksi)*(ksi+2*eta),   0.25*(-1+ksi)*(ksi-2*eta),...
                      0.5*(-1+ksi)*(1+ksi),                -(1+ksi)*eta,       -0.5*(-1+ksi)*(1+ksi),               (-1+ksi)*eta ];
        end

    case {'Q4', 'AHMAD4'}
        dN = zeros(2,4,ngauss);
        for igauss = 1:ngauss
            
            ksi = pointArray(igauss,1);
            eta = pointArray(igauss,2);

            dN(:,:,igauss) = [  % derivadas respecto de ksi
                -0.25*(1 - eta),  0.25*(1 - eta), 0.25*(1 + eta), -0.25*(1 + eta)
                % derivadas respecto de eta
                -0.25*(1 - ksi), -0.25*(1 + ksi), 0.25*(1 + ksi),  0.25*(1 - ksi) ];
        end
        
    case 'H8' % Segun Cook pag 217.
        dN = zeros(3,8,ngauss);
        for igauss = 1:ngauss
            
            ksi  = pointArray(igauss,1);
            eta  = pointArray(igauss,2);
            zeta = pointArray(igauss,3);
            
            dN(:,:,igauss) = [ % derivadas respecto de ksi
                eta/8 - zeta/8 + (eta*zeta)/8 - 1/8,  eta/8 + zeta/8 - (eta*zeta)/8 - 1/8,...
                zeta/8 - eta/8 + (eta*zeta)/8 - 1/8,  - eta/8 - zeta/8 - (eta*zeta)/8 - 1/8,...
                zeta/8 - eta/8 - (eta*zeta)/8 + 1/8, (eta*zeta)/8 - zeta/8 - eta/8 + 1/8,...
                eta/8 - zeta/8 - (eta*zeta)/8 + 1/8,  eta/8 + zeta/8 + (eta*zeta)/8 + 1/8 
                % derivadas respecto de eta
                ksi/8 - zeta/8 + (ksi*zeta)/8 - 1/8, ksi/8 + zeta/8 - (ksi*zeta)/8 - 1/8,...
                (ksi*zeta)/8 - zeta/8 - ksi/8 + 1/8, zeta/8 - ksi/8 - (ksi*zeta)/8 + 1/8,...
                - ksi/8 - zeta/8 - (ksi*zeta)/8 - 1/8, zeta/8 - ksi/8 + (ksi*zeta)/8 - 1/8,...
                ksi/8 - zeta/8 - (ksi*zeta)/8 + 1/8, ksi/8 + zeta/8 + (ksi*zeta)/8 + 1/8
                % derivadas respecto de zeta
                (eta*ksi)/8 - ksi/8 - eta/8 + 1/8, eta/8 + ksi/8 - (eta*ksi)/8 - 1/8,...
                ksi/8 - eta/8 + (eta*ksi)/8 - 1/8, eta/8 - ksi/8 - (eta*ksi)/8 + 1/8,...
                ksi/8 - eta/8 - (eta*ksi)/8 + 1/8, eta/8 - ksi/8 + (eta*ksi)/8 - 1/8,...
                - eta/8 - ksi/8 - (eta*ksi)/8 - 1/8,  eta/8 + ksi/8 + (eta*ksi)/8 + 1/8];
        end
        
      case 'H8ADINA' % Segun manual de ADINA.
        dN = zeros(3,8,ngauss);
        for igauss = 1:ngauss
            
            ksi  = pointArray(igauss,1);
            eta  = pointArray(igauss,2);
            zeta = pointArray(igauss,3);
            
            dN(:,:,igauss) = [ % derivadas respecto de ksi
                -((eta - 1)*(zeta + 1))/8,  ((eta - 1)*(zeta + 1))/8,...
                -((eta - 1)*(zeta - 1))/8,  ((eta - 1)*(zeta - 1))/8,...
                 ((eta + 1)*(zeta + 1))/8, -((eta + 1)*(zeta + 1))/8,...
                 ((eta + 1)*(zeta - 1))/8, -((eta + 1)*(zeta - 1))/8 
                % derivadas respecto de eta
                -(ksi/8 + 1/8)*(zeta + 1),  (ksi/8 - 1/8)*(zeta + 1),...
                -(ksi/8 - 1/8)*(zeta - 1),  (ksi/8 + 1/8)*(zeta - 1),...
                 (ksi/8 + 1/8)*(zeta + 1), -(ksi/8 - 1/8)*(zeta + 1),...
                 (ksi/8 - 1/8)*(zeta - 1), -(ksi/8 + 1/8)*(zeta - 1)
                % derivadas respecto de zeta
                 -(ksi/8 + 1/8)*(eta - 1),  (ksi/8 - 1/8)*(eta - 1),...
                 -(ksi/8 - 1/8)*(eta - 1),  (ksi/8 + 1/8)*(eta - 1),...
                  (ksi/8 + 1/8)*(eta + 1), -(ksi/8 - 1/8)*(eta + 1),...
                  (ksi/8 - 1/8)*(eta + 1), -(ksi/8 + 1/8)*(eta + 1)];
        end
        
        case 'H8FACU' % Segun criterio Elements_sorter.
        dN = zeros(3,8,ngauss);
        for igauss = 1:ngauss
            
            ksi  = pointArray(igauss,1);
            eta  = pointArray(igauss,2);
            zeta = pointArray(igauss,3);
            
            dN(:,:,igauss) = [% derivadas respecto de ksi
               -((eta + 1)*(zeta + 1))/8,  ((eta - 1)*(zeta + 1))/8,...
               -((eta - 1)*(zeta + 1))/8,  ((eta + 1)*(zeta + 1))/8,...
                ((eta + 1)*(zeta - 1))/8, -((eta - 1)*(zeta - 1))/8,...
                ((eta - 1)*(zeta - 1))/8, -((eta + 1)*(zeta - 1))/8  
                             % derivadas respecto de eta
               -(ksi/8 - 1/8)*(zeta + 1),  (ksi/8 - 1/8)*(zeta + 1),...  
               -(ksi/8 + 1/8)*(zeta + 1),  (ksi/8 + 1/8)*(zeta + 1),...
                (ksi/8 - 1/8)*(zeta - 1), -(ksi/8 - 1/8)*(zeta - 1),...
                (ksi/8 + 1/8)*(zeta - 1), -(ksi/8 + 1/8)*(zeta - 1) 
                            % derivadas respecto de zeta
               -(ksi/8 - 1/8)*(eta + 1),  (ksi/8 - 1/8)*(eta - 1),...
               -(ksi/8 + 1/8)*(eta - 1),  (ksi/8 + 1/8)*(eta + 1),...
                (ksi/8 - 1/8)*(eta + 1), -(ksi/8 - 1/8)*(eta - 1),...
                (ksi/8 + 1/8)*(eta - 1), -(ksi/8 + 1/8)*(eta + 1)];
        end
        
end

