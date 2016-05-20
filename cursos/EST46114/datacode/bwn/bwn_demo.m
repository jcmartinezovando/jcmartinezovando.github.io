% DEMO

clear;
echo off;

addpath = 'C:\JCMO.Trabajo\@Mis.Cursos\2016-I_Inferencia Bayesiana en Alta Dimension\_data&code\bwn';
cd('C:\JCMO.Trabajo\@Mis.Cursos\2016-I_Inferencia Bayesiana en Alta Dimension\_data&code\bwn');

phis=[0.02
    0.02969
    -0.68
    0.0068
    0.305 ];
raices=roots(phis)

% ======================================================================
% Simulacion del Proceso
p=max(size(phis));
sigma=5;
Y0=zeros(p,1);
T=160;
% Generar el proceso
[Y] = ar_sim(p,phis,sigma,T,Y0);

%figure(1)
%plot(Y);
T=length(Y);
N=15;  % ultimos valores para predecir
Y_fut=Y(T-N+1:1:T);
Y_pred_mcmc=cell(N,1);       % repositorio de las muestras de los valores predictivos 
Y_pred_media=zeros(N,1);

cont=0;
t = 1;
for t=T-N:T-1
    fprintf('           ##                      -------------------------                   Estamos en el tiempo     %d   de        %  .\n',t,T);
    cont=cont+1;
    Ytemp=Y(1:1:t);
    [mcmc,y_prediccion,y_pred_media]=bwn(Ytemp);
    Y_pred_mcmc{cont}=y_prediccion;
    Y_pred_media(cont)=y_pred_media;
    clear mcmc;    clear y_prediccion; clear y_pred_media;     % LIBERAMOS espacio necesario
end

figue(10)
t=(T-N+1:1:T);
plot(t,Y_fut,Y_pred_media,'r-','b:');
hold off;
