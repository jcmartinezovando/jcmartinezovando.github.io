function [Ypred] = bwn_pred(Y,mcmc,op)
%   Calcula una muestra de la deistribuci´on predictiva un paso adelante usando el modelo de onduletas radiales
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

M=op.mcmc;
Ypred=zeros(M,1);
T=length(Y);
yf=Y(T:-1:T-op.p+1)';
for m=1:M
    Xwnf=bwnmat(yf,mcmc.ks{m},mcmc.cent{m},mcmc.dilat{m},op.wavelet);
    Xwnf=[1 yf Xwnf];
    Ypred(m)=mvnrnd(Xwnf*mcmc.omegas{m}/200,sqrt(mcmc.taus{m}.^-1),1);
end
%
%   --  FIN --