function Meshplot(elementos,nodos,color,verbatim)
% MESHPLOT  Graficador de mallas ADAPTADO PARA 3D
% elementos: matriz de conectividades.
% nodos:     matriz de coordenadas nodales.
% color:     string que especifica el color para los bordes de los elementos.

sups = reshape(elementos(:,[1 2 3 4  5 6 7 8  3 4 8 7  2 1 5 6  1 5 8 4  2 6 7 3])',4,[])'; 

h1 = patch('Faces',sups,'Vertices',nodos);
set(h1,'EdgeColor',color,'FaceColor','none');

set(gca,'XTick',[],'YTick',[],'ZTick',[],'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0])
daspect([10 10 10])
hold on
if verbatim
    for n=1:size(nodos,1)        
        text(nodos(n,1),nodos(n,2),nodos(n,3),num2str(n),'FontSize',12,'Color','r')
    end
    
    for e=1:size(elementos,1)
        text(mean(nodos(elementos(e,:),1)),mean(nodos(elementos(e,:),2)),mean(nodos(elementos(e,:),3)),num2str(e),'FontSize',12,'Color','b')
    end
end
xlabel('X direction','Color','k');
ylabel('Y direction','Color','k');
zlabel('Z direction','Color','k');
end
