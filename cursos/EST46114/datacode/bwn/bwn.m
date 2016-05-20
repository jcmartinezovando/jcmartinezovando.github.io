function [mcmc,y_prediccion,y_pred_media,op] = bwn(Y,op)
%       Estima el modelo de regresión semiparametrica basado en onduletas radiales
%       
%
%       ENTRADAS:
%                  Y   -   vector de datos de la serie de tiempo de longitud ' n '
%                 op   -   opciones del modelo en forma de  ' lista '
%                           op.lin          -       variable indicadora:  
%                                                               0 - no incluye termino lienal, 
%                                                               1 - incluye termino lineal
%                           op.kmax     -       numero maximo de funciones rediales
%                           op.p        -       orden de autorregresion del proceso
%                           op.mcmcm    -   tamaño de muestras despu´es del periodo de calentamiento
%                           op.calentamiento      -       longitud del periodo de calentamiento
%                           op.precision    -   parametro de precision de los coeficientes de regression
%                           op.gamma1       -       parametro de la distribucion inicial de los grados de libertad
%                           op.gamma2       -       segundo parametro  "                    "                   "
%                                                               gl  ~ Gamma(gamma1,gamma2)    
%                           op.alfa1    -   parametro de la varianza de los errores
%                           op.alfa2    -   parametro de la varianza de los errores
%                                                               tau=1/sigma^ 2  ~  Gamma(alfa1,alfa2)
%                            op.prob    -   vector con cuatro entradas de las probabilidades de transicion de cuatro movimientos
%                                                   op.prob=(a,b,r,t)
%                                                                   a    -  prob. de proponer adicionar una base
%                                                                   b    -  prob. de proponer eliminar una base
%                                                                    r    -  prob. de alterar el parametro de escala
%                                                                    t    -  prob. de cambiar el parametro de traslacion  
%                           op.wavelet      -       tipo de onduelta que usaremos
%                           op-graficas     -       indicadora de graficas o no
%                           op.estandar     -       opcion de estandarizar la matriz de covariables
%
%       SALIDAS:
%                           mcmc        -       lista de variables generadas por el muesterador
%
%   Compatible con MATLAB  5.0
%

%   --------------------------------------                          OPCIONES
%
if nargin == 2
else % usamos las opciones estandar
    op.lin=1;
    op.kmax=50;
    op.kmin=33;	
    op.p=3;
    op.mcmc=10000;
    op.calentamiento=20000;
    op.precision=0.1;
    op.alfa1=10^-2;
    op.alfa2=10^ -2;
    op.gamma1=10^-3;
    op.gamma2=10^-3;
    op.prob=[0.15,0.15,0.35,0.35];
    op.wavelet='marr';
    op.K=3;
    op.Q=1.1;
    op.graficas=1;
    op.estandar=1;
end
%   -------------------------------------                       Preparamos los datos
% Arreglo de la serie en forma matricial
p=op.p;
n=length(Y); 
K=op.K; Q=op.Q;
y=Y(n:-1:p+1);             
x=flipud(hankel(Y(n-1:-1:p),Y(p:-1:1))); 
y=Y(p+1:1:n);
n=length(y);
xf=y(n:-1:n-p+1)';

% clear Y;            % liberamos espacio
[n d]=size(x);
yty=y'*y;

if op.estandar
   for i=1:d
      mx(i) = mean(x(:,i));
      sx(i) = std(x(:,i));
      x(:,i) = (x(:,i)-mx(i)) / sx(i);
      xf(i)=(xf(i)-mx(i)) / sx(i);
  end
else
   mx=ones(1,p); sx=ones(1,p);
end
mcmc.mx=mx;
mcmc.sx=sx;

%   -------------------------------------                       inicializamos las variables por guardar
mcmc.ypreds=cell(op.mcmc,1);                % muestra de los valores predictivos
mcmc.yajust=cell(op.mcmc,1);                % muestra de los valores ajustados
mcmc.Xwns=cell(op.mcmc,1);                  % matriz de diseño
mcmc.omegas=cell(op.mcmc,1);                % muestra de los coeficientes de regresion
mcmc.omegas_media=cell(op.mcmc,1);      % muestra de las medias de los coef. de  regresion
mcmc.taus=cell(op.mcmc,1);                      % muestra de tau
mcmc.ks=cell(op.mcmc,1);                            % muestra del numero de bases
mcmc.tSs=cell(op.mcmc,1);                           % muestra de los grados de liebrtad del modelo
mcmc.cent=cell(op.mcmc,1);
mcmc.dilat=cell(op.mcmc,1);
% registro de los moviemienos
mcmc.propuestas=zeros(4,1);                         
mcmc.aceptaciones=zeros(4,1);

% determinamos los parametros de dilatacion y traslacion iniciales
D=K.*ones(op.kmax,1);                                               % dilatacion
dd=length(D);
for l=1:dd
    D(l)=D(l)^1;
end
Mu=zeros(d,op.kmax);                                                  % traslacion
            % seleccionamos 'op.kmax' vectores de la matriz x (sin remplazo) 
[mue] = muest_sr((1:1:n),op.kmax);
for j=1:op.kmax
    Mu(:,j)=Q.*x(mue(j),:)';
end

% Inicializamos la matriz de diseño de la Base de Onduletas
% en este caso inicializaremos la cadena en la mitad de las bases de oinduletas maxima 
if op.lin == 1      % incluimos un componente lineal   
    kmin=op.kmin;
    kmax=op.kmax;
    % el numero incial del numero de bases de ondueltas
    k=kmax;      % numero incial de bases
    Xwn=zeros(n,d+1+k);
    Xwn(:,1:d+1)=[ones(n,1) x]; % esta pendiente el calculo de las demas entradas de la matriz
    Xwn(:,d+2:d+1+k)=bwnmat(x,k,Mu,D,op.wavelet);
else
    kmin=op.kmin;
    kmax=op.kmax;
    % el numero incial del numero de bases de ondueltas
    k=kmax;     % numero incial de bases
    Xwn=zeros(n,k+1);
    Xwn(:,1)=ones(n,1); % esta pendiente el calculo de las demas entradas de la matriz
    Xwn(:,2:k+1)=bwnmat(x,k,Mu,D,op.wavelet);
end
%  ------------------------------------------------------------------------------------------------------------------------------------
%        RJMCMC             -       inicio
%                   inicializamos algunas de las variables
% probabilidades de transicion de modelos
pn=op.prob(1); pm=op.prob(1)+op.prob(2); pd=op.prob(1)+op.prob(2)+op.prob(3);
 % variables de conteo
conteo=0;                 %   contador para el numero de movientos de la cadena 
muestra=0;              % contador de muestras mcmc recopiladas
k_prop=k;                %   inicializamos la indicadora del orden del modelo ' propuesto '
% registro de propuestas y aceptacion
propuesta=zeros(4,1);
aceptacion=zeros(4,1);
indicadora=0;       % variable indicadora auxiliar para el registro de propuestas y aceptaciones
%   ----------------------------            Tomamos una muestra de tau, omega; y calculamos log-verosimil 
%                                                     del modelo inicial
tau=(op.alfa2^-1)*randgamma_mat(op.alfa1,1,1);
[log_mv,omega,omega_media,sumcuad,tS,tS_prob] = bwnmv(y,yty,Xwn(:,1:1+d+k),tau,...
   op.precision,op.alfa1,op.alfa2,...
   op.gamma1,op.gamma2);
fprintf('Iniciamos MCCM - SR...\n');
%       --------------------------              Desarrollo de la cadena
while muestra < op.mcmc
    conteo=conteo+1;
    % indicadoras de los movimientos
    nacimiento=0;       % indicadora si aumenamos en una base
    muerte=0;               % indicadora si eliminamos una base
    dilatacion=0;           % indicadora de la dilatacion
    traslacion=0;           % indicadora de la traslacion
    %       -------------------------           Hacemos una copia del modelo actual en las variables 'prop'
    %                                            si rechazamos el movimiento, estas seran las variables propuestas
    omega_prop=omega;
    Xwn_prop=Xwn(:,1:k+d+1);
    omega_temp=omega_prop;
    Xwn_temp=Xwn_prop;
    k_temp=k;
    %       -------------------------           Determinamos el tipo de movimiento propuesto
    u1=rand;
    if 0 < u1 <= pn  % añadimos una base de onduletas
        nacimiento=1;
        indicadora=1;
        if k == kmax
            nacimiento=0;
            traslacion=1;
            indicadora=4;
        end
    elseif pn < u1 <= pm    % eliminamos una base de onduletas
        muerte=1;
        indicadora=2;
        if k == kmin    % eliminamos (aleatoriamente) una base de onduletas
            muerte=0;
            dilatacion=1;
            indicadora=3;
        end
    elseif pm < u1 <= pd     %   cambiamos de escala una de las bases
        dilatacion=1;
        indicadora=3;
    else        %  cambiamos alguno de los parametros de localizacion
        traslacion=1;
        indicadora=4;
    end    
    % guardamos el modelo que proponemos
    propuesta(indicadora)=propuesta(indicadora)+1;    
    % 
    % --------------     Actualizacion del modelo propuesto de acuerdo al tipo de movimiento propuesto
    %    
    if nacimiento       %           -----------------------         NACIMIENTO
        k_prop=k+1;
        % tenemos que actualizar la matriz de diseño de acuerdo al modelo porpuesto (aleatoria)
        Xwn_prop(:,1+d+k_prop)= bwnnac(x,k_prop,Mu,D,op.wavelet);
    elseif muerte           %       -----------------------             MUERTE
        k_prop=k-1;
        % elegimos aleatoriamente una de las 'k' bases para eliminarla
        k_elim=ceil(rand*k);
        Xwn_prop(:,1+d+k_elim)=[];
    elseif dilatacion       %       -----------------------                 DILATACION
        k_prop=k;
        % elegimos aleatoriamente una de la 'k' bases para modificar su 'dilatacion'
        k_dil=ceil(rand*k);
        [Xwn_prop(:,k_dil), dil_prop]= bwndil(x,k_dil,Mu,D,op.wavelet);
    elseif traslacion       %       ------------------------                TRASLACION
        k_prop=k;
        % elegimos aleatoriamente una de la 'k' bases para modificar su 'traslacion'
        k_tras=ceil(rand*k);
        [Xwn_prop(:,k_tras+d+1), Mu_prop]= bwntras(x,k_tras,Mu,D,op.wavelet);
        % Mu_prop es un vector de dimension 'dx1'
    end
    %
    %   ----------------------         Calculo de la log-verosimilitud (funcion auxiliar) del modelo propuesto 
    %
    [log_mv_prop,omega_prop,omega_media_prop,sumcuad_prop,tS_prop,tS_prob_prop] = ...
                                                                 bwnmv(y,yty,Xwn_prop(:,1:1+d+k_prop),tau,op.precision,...
                                                                                op.alfa1,op.alfa2,op.gamma1,op.gamma2);
    %
    %   -----------------------                 Aceptacion o cancelacion del modelo 
    %
    if rand < exp(log_mv_prop-log_mv)*(tS_prob_prop/tS_prob)
        %   --  ACEPTAMOS   --  el movimiento
        omega=omega_prop;
        omega_media=omega_media_prop;
        sumcuad=sumcuad_prop;
        k=k_prop;
        Xwn(:,1:1+d+k)=Xwn_prop;
        aceptacion(indicadora)=aceptacion(indicadora)+1;
        log_mv=log_mv_prop;
        tS_prob=tS_prob_prop;
        if dilatacion
            D(k_dil)=dil_prop;
        elseif traslacion
            Mu(:,k_tras)=Mu_prop;
        end
    else
        
    end
    %
    %   -----------------------                 Obtenemos un valor de 'tau' (precision) y un valor predictivo un paso adelante
    %
    tau=(1./(op.alfa2+sumcuad)/2)*randgamma_mat((op.alfa1+n)/2,1,1);
    Xwn_pred=bwnmat(xf,k,Mu,D,op.wavelet);
    Xwn_pred=[1 xf Xwn_pred];
    y_pred=mvnrnd(Xwn_pred*omega,sqrt(tau^-1),1);
    y_ajust=Xwn(:,1:1+d+k)*omega;
    %
    %   -----------------------                 Verificamos el periodo de calentamiento
    %
    if conteo > op.calentamiento
        % Empezamos a guardar las muestras en 'mcmc'
        muestra=muestra+1;
        mcmc.ypreds{muestra}=y_pred;
        mcmc.yajust{muestra}=y_ajust;
        mcmc.omegas{muestra}=omega;
        mcmc.taus{muestra}=tau;
        mcmc.ks{muestra}=k;
        mcmc.tSs{muestra}=tS;
        mcmc.Xwns{muestra}=Xwn(:,1:1+d+k);
        mcmc.cent{muestra}=Mu;
        mcmc.dilat{muestra}=D;
    end
    % Imprimimos el conteo de iteraciones
    if rem(conteo,10) == 0
        if muestra == 0
            fprintf('Calentamiento       %d  /  %d  \n',conteo,op.calentamiento);
        else
            fprintf('Iteracion       %d  /  %d  \n',muestra,op.mcmc);
        end
    end
    %
end    % fin del while
mcmc.propuestas=propuesta;
mcmc.aceptaciones=aceptacion;
%
y_prediccion=zeros(op.mcmc,1);
y_pred_erg=zeros(op.mcmc,1);
t=(1:1:op.mcmc);
for i=1:op.mcmc
    y_prediccion(i)=mcmc.ypreds{i};
end
for i=1:op.mcmc
    y_pred_erg(i)=([ones(1,i) zeros(1,op.mcmc-i)]*y_prediccion)/i;
end
%
%           -----------------------------               Graficas
%
if op.graficas
    %   --      grafica de la evolucion del modelo
    modelo=zeros(op.mcmc,1);
    t=(1:1:op.mcmc);
    for i=1:op.mcmc
        modelo(i)=mcmc.ks{i};
    end
    figure(1)
    subplot(311)
    hold on;
    plot(t,modelo,'r');
    axis([0 op.mcmc 0 op.kmax+1]);
    ylabel('Modelo','fontsize',10); 
    xlabel('ite','fontsize',10)    
    subplot(312);
    hist(modelo,1:1:op.kmax)
    ylabel('p(M_{k})','fontsize',10);
    xlabel('Modelo','fontsize',10);
    hold on;
    subplot(313)
    autocorr(modelo);
    grid off;
    xlabel('Retraso');
    ylabel('FAC muestral');
    title('');
    hold off;
    clear modelo;
        %   --      grafica de la evolucion de los valores predictivos 
   
    figure(2)
    subplot(221)
    hold on;
    plot(t,y_prediccion,'r');
    axis([0 op.mcmc min(y_prediccion)-0.5 max(y_prediccion)+0.5]);
    ylabel('y_{futura}','fontsize',10); 
    xlabel('ite','fontsize',10)    
    subplot(222)
    hold on;
    plot(t,y_pred_erg,'r');
    axis([0 op.mcmc min(y_prediccion)-0.5 max(y_prediccion)+0.5]);
    ylabel('Promedio Ergodico','fontsize',10); 
    xlabel('ite','fontsize',10)    
    subplot(223);
    hist(y_prediccion,100)
    ylabel('p(y_{futura})','fontsize',10);
    xlabel('y_{futura}','fontsize',10);
    y_pred_media=mean(y_prediccion);
    hold on;
    subplot(224)
    autocorr(y_prediccion);
    grid off;
    xlabel('Retraso');
    ylabel('FAC muestral');
    title('');
    hold off;
    clear y_pred_erg;
    hold off;
else
end
%
%       ---         FIN     ---