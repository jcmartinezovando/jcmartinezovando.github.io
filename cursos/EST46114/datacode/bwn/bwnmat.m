function [Wf] = bwnmat(X,k,Mu,D,wavelet)
% Construye una matriz de diseño de bases de onduletas con la funcion 'wavelet'
%
%       NOTA:     Esta funcion es invocada dentro de 'bwn.m'
%
%       ENTRADAS:
%               X       -       matriz de variables regresoras de dimension 'nxp'
%               k       -       numero de funciones base de onduletas por construir
%               Mu      -       matriz de parametros de traslacion de 'pxk'
%               D       -       vectorde parametros de dilatacion de 'kx1'
%           wavelet     -       tipo (familia) de onduleta
%
%       SALIDAS:
%               Wf       -       matriz de diseño de las funciones de onduletas
%

[n,p]=size(X);
W=zeros(n,k);
Wf=zeros(n,k);
% construimos la matriz de traslacion y dilatacion
for j=1:k
    for i=1:n
        w_temp=X(i,:)-Mu(:,j)';
        norm2=(w_temp)*(w_temp)';
        W(i,j)=D(j)*(norm2^0.5);
    end
end
% evaluamos la funcion 'wavelet'
switch lower(wavelet)
case 'marr'
    Wf=(p-(W.^2)).*exp(-0.5*(W.^2));
end

%      --   FIN --