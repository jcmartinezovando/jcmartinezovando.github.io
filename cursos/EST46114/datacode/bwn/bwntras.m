function [Wf , mu]= bwntras(X,k_tras,Mu,D,wavelet)
% Traslada aleatoriamenet una columna de Xwn  con base en X (original)
%
%       NOTA:     Esta funcion es invocada dentro de 'bwn.m'
%
%       ENTRADAS:
%               X       -       matriz de variables regresoras de dimension 'nxp'
%               k_tras      -      base de Xwn que se trasladara
%               Mu      -       matriz de parametros de traslacion de 'pxk'
%               D       -       vectorde parametros de dilatacion de 'kx1'
%           wavelet     -       tipo (familia) de onduleta
%
%       SALIDAS:
%               Wf       -       matriz (columna) de diseño de las funciones de onduletas
%

[n,p]=size(X);
W=zeros(n,1);
Wf=zeros(n,1);
% elegimos aleatoriamente un nuevo centroide de la base
ind=ceil(rand*n);
mu=X(ind,:)';
% construimos la matriz de traslacion y dilatacion
for i=1:n
    w_temp=X(i,:)-mu';
    norm2=(w_temp)*(w_temp)';
     W(i,1)=D(k_tras)*(norm2^0.5);
end
% evaluamos la funcion 'wavelet'
switch lower(wavelet)
case 'marr'
    Wf=(p-(W.^2)).*exp(-0.5*(W.^2));
end

%      --   FIN --