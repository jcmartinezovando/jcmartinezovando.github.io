function [Wnac] = bwnnac(X,k_prop,Mu,D,wavelet)
% Genera una columna de base de onduletas indice 'k_prop'
%
%       NOTA:     Esta funcion es invocada dentro de 'bwn.m'
%
%       ENTRADAS:
%               X       -       matriz de variables regresoras de dimension 'nxp'
%               k_prop      -      base de Xwn que se trasladara
%               Mu      -       matriz de parametros de traslacion de 'pxk'
%               D       -       vectorde parametros de dilatacion de 'kx1'
%           wavelet     -       tipo (familia) de onduleta
%
%       SALIDAS:
%               Wf       -       matriz (columna) de diseño de las funciones de onduletas
%

[n,p]=size(X);
W=zeros(n,1);
Wnac=zeros(n,1);
% construimos la matriz de traslacion y dilatacion
for i=1:n
    w_temp=X(i,:)-Mu(:,k_prop)';
    norm2=(w_temp)*(w_temp)';
     W(i,1)=D(k_prop)*(norm2^0.5);
end
% evaluamos la funcion 'wavelet'
switch lower(wavelet)
case 'marr'
    Wnac=(p-(W.^2)).*exp(-0.5*(W.^2));
end

%      --   FIN --