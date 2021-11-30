clear all
close all
clc
%% Codigo para analizar evolucion de T de resolución Vs Cantidad de superelmentos que se usan en el modelo:
tRes = [0.01 0.012 0.02 0.031 0.122 0.054 0.144 0.066 0.086 0.081 0.166 0.182 0.141 0.181 0.272 0.306 0.307 0.381 0.362 0.355 0.467];
nDof = [600  744   888  1032  1176  1248  1320  1392  1464  1536  1608  1752  1896  2040  2184  2328  2472  2616  2760  2904  3048];

plot(nDof,tRes);
grid on
title('Tiempos de resolución vs N° de Dofs');
xlabel('N° de dofs Totales');
ylabel('Tiempo de resolución en [s]');
