function [log_mv,omega,omega_media,sumcuad,tS,tS_prob] = bwnmv(y,yty,X,tau,precision,alfa,beta,gamma1,gamma2)
%  Calcula la pseudo log-maxima verosimilitud del modelo necesaria para el Factor de Bayes
%  en la evaluacion de la probabilidad de aceptacion del movimiento
%
%  NOTA:    La salida de esta funcion es util para calcular la probabilidad de aceptacion del movi-
%                   miento propuesto en RJMCMC.
%
%   ENTRADAS:
%               y       -       vector de respuestas
%               yty     -       suma de cuadrados de y
%               X       -       matriz de diseño del modelo
%               tau     -       muestra de la precision de los errores
%               precision    -   precision sobre 'omega'
%               alfa, beta    -   parametros de la distribucion sobre 'tau'
%
%   SALIDAS:
%               log_mv      -       pseudo log- Verosimilitud del modelo
%               omega       -       muestra de la distribucion final de omega
%               omega_media      -   media de la distribucion final de 'omega'  
%               sumcuad     -   suma de cuadrados final
%               tS     -       trasa de la matriz S (suavizamiento)
%               tS_prob     -       pesudo-probabilidad inicial sobre los grados de libertad
%

[n p]=size(X);
Pini=precision*eye(p);
Pfin=X'*X+Pini;
Vfin=inv(Pfin);
omega_media=Vfin*X'*y;
T=chol(Vfin);               % matriz trianguilar superior en la descomposicion Cholesky
% suma de cuadrados final
sumcuad=yty-omega_media'*Pfin*omega_media;
% muestra de la distribucion final de 'omega'
V=(tau^ -1).*Vfin;
omega=mvnrnd(omega_media,V,1)';
% calculamos la log-verosimilitud del modelo
log_detfin=sum(log(diag(T)));
log_detini=-0.5*p*log(precision);
% usaremos la suma de cuadrados que calculamos previamente
% util para calcular el factor de Bayes del modelo propuesto
log_mv=log_detfin-log_detini-(0.5*(n+beta))*log(0.5*(alfa+sumcuad));
% calculamos la matriz de suavizamiento y su traza
S=X*Vfin*X';
tS=trace(S);
% calculamos la pesudo-probabilidad inicial (salvo por una constante), 
% util para clacular los momios inciales del modeo propuesto
tS_prob=tS^(gamma1-1)*exp(-gamma2*tS);
%
%    --     FIN     --

