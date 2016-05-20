function [Y] = ar_sim(p,phi,sigma,T,Y0)
% Simula un proceso autorregresivo est'atico de orden p
% Entradas:
%          p   		- orden del proceso
%          phi    	- vector de coeficientes de autorregresion de px1
%          sigma   	- varianza del proceso
%          T   		- longitud del proceso
%          Y0  		- vector de trayectorias iniciales del proceso de imensin px1 
%
% Salidas: 
%          Y  		- vector de px1 del proceso simulado
%
if nargin < 5 
    error('No hay suficientes argumentos.');
end

if  max(size(phi)) ~= p
    error('El vector de coeficientes no coincide con el orden del proceso. ')
end

if  max(size(Y0)) ~= p
    error('La dimesnion de la trayectoria inicial no coincide con el orden del proceso. ')
end

if  sigma <= 0
    error('La varianza no esta especificada correctamente. ')
end

Yt=zeros(p,1);                              % es el vector de coeficientes de autorregresion temporal
Y=zeros(T,1);                               
Tt=Y0;

for t = 1:T
    Y(t)=Yt'*phi+normrnd(0,sqrt(sigma));
    if t <= p
        Yt=[Y(t:-1:1);Yt(p:-1:t+1)];
    else
        Yt=Y(t:-1:t-p+1);        
    end
end
