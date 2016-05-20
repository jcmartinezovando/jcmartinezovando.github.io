function [verosimil] = bwn_verint(Y,mcmc,op)
%      Calcula una muestra de la distribución predictiva un paso adelante usando el modelo de onduletas radiales
%       
%   NOTA: Requiere las salidas "mcmc" de la funcion 'bwn.m'
%
%       ENTRADA
%                  Y   -   vector de datos de la serie de tiempo de longitud ' n '
%          mcmc        -       lista de variables generadas por el muesterador
% 
%       SALIDA
%               Ypred   -   muestra de la distribucion predictiva un paso adelante
%

%   -------------------------------------                       Preparamos los datos
% Arreglo de la serie en forma matricial
p=op.p;
T=length(Y);
yy=Y(p+1:T);
N=length(yy);
suma=0;
%temp2=1;
for m=1:op.mcmc
    if rem(m,20)==0
        fprintf('Estamos en la iteracion %d de %d.\n',m,op.mcmc);
    end
    prod=1;
    for k=1:N
        prod=prod*normpdf(yy(k), mcmc.Xwns{m}(k,:)*mcmc.omegas{m},sqrt(mcmc.taus{m}^-1));
    end
    suma=suma+inv(prod);
end
suma=suma/op.mcmc;
verosimil=inv(suma);
clear suma; clear prod;
%
%   --  FIN --
