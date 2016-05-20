function [Wf , dil] = bwndil(X,k_dil,Mu,D,wavelet)
% Dilata aleatoriamenet una columna de Xwn  con base en X (original)
%
%       NOTA:     Esta funcion es invocada dentro de 'bwn.m'
%
%       ENTRADAS:
%               X       -       matriz de variables regresoras de dimension 'nxp'
%               k_dil      -     columna de Xwn que se dilatara
%               Mu      -       matriz de parametros de traslacion de 'pxk'
%               D       -       vectorde parametros de dilatacion de 'kx1'
%           wavelet     -       tipo (familia) de onduleta
%
%       SALIDAS:
%               Wf       -       matriz (columna) de diseño de las funciones de onduletas
%               dil         -       dilatacion propuesta
%

[n,p]=size(X);
W=zeros(n,1);
Wf=zeros(n,1);
dil=0;
dil=exp(log(D(k_dil))+mvnrnd(0,1.5,1));
% construimos la matriz de traslacion y dilatacion
for i=1:n
    w_temp=X(i,:)-Mu(:,k_dil)';
    norm2=(w_temp)*(w_temp)';
     W(i,1)=dil*(norm2^0.5);
end
% evaluamos la funcion 'wavelet'
switch lower(wavelet)
case 'marr'
    Wf=(p-(W.^2)).*exp(-0.5*(W.^2));
end
%
%      --   FIN --