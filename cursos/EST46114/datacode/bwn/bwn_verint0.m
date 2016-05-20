function [MC] = bwn_verint(Y,op,M)
%     Esta funcion calcula el estimador de Monte Carlo de la verosilitud integrada del modelo de bases radiales de ondueltas.
%
%     NOTA:   'op' debe ser la misma que se utiliza en 'bwm.m'.
%       
%
%     ENTRADAS:
%                  Y   -   vector de datos de la serie de tiempo de longitud ' n '
%                 M  -  tamaño de la muestra del estimador de Monte Carlo
%               op   -   opciones del modelo en forma de  ' lista '
%                           op.lin          -       variable indicadora:  
%                                                               0 - no incluye termino lienal, 
%                                                               1 - incluye termino lineal
%                           op.kmax     -       numero maximo de funciones rediales
%                           op.kmin     -       numero minimo de funciones rediales
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
%
%       SALIDAS:
%                           MC        -       Estimador de Monte Carlo de la veroslimiltud integrada
%

%   -------------------------------------                       Preparamos los datos
% Arreglo de la serie en forma matricial
p=op.p;
n=length(Y); 
K=op.K; Q=op.Q;
y=Y(n:-1:p+1);             
x=flipud(hankel(Y(n-1:-1:p),Y(p:-1:1))); 
y=Y(p+1:1:n);
n=length(y);
xf=y(n-p+1:1:n)';
K=op.K;
% clear Y;            % liberamos espacio
[n d]=size(x);
yty=y'*y;
H=length(y);

% generamos una muestra de tamaño M de los modelos
m= rndud(op.kmin,op.kmax,M);

% determinamos los parametros de dilatacion y traslacion iniciales
D=K.*ones(op.kmax,1);                                               % dilatacion
dd=length(D);
for l=1:dd; D(l)=D(l)^l; end;
Mu=zeros(d,op.kmax);                                                  % traslacion
            % seleccionamos 'op.kmax' vectores de la matriz x (sin remplazo) 
[mue] = muest_sr((1:1:n),op.kmax);
for j=1:op.kmax
    Mu(:,j)=Q.*x(mue(j),:)';
end

% Inicializamos la matriz de diseño de la Base de Onduletas en el numero maximo de onduletas
if op.lin == 1      % incluimos un componente lineal   
    % el numero incial del numero de bases de ondueltas
    k=op.kmax;      % numero incial de bases
    Xwn=zeros(n,d+1+k);
    Xwn(:,1:d+1)=[ones(n,1) x]; % esta pendiente el calculo de las demas entradas de la matriz
    Xwn(:,d+2:d+1+k)=bwnmat(x,k,Mu,D,op.wavelet);
else
    % el numero incial del numero de bases de ondueltas
    k=op.kmax;     % numero incial de bases
    Xwn=zeros(n,k+1);
    Xwn(:,1)=ones(n,1); % esta pendiente el calculo de las demas entradas de la matriz
    Xwn(:,2:k+1)=bwnmat(x,k,Mu,D,op.wavelet);
end

% generamos una muestra de los coeficientes de regresion
betam=cell(M,1);
gamam=NaN.*zeros(M,1);
for i = 1:M
    betam{i,1}=mvnrnd(zeros(m(i)+d+1,1),inv(op.precision.*eye(1+d+m(i))),1)';
    gamam(i)=(op.alfa2).*randgamma_mat(op.alfa1,1,1);
end

SUMA=0;
for i=1:M
    constante1=(2*pi)^(1+d+m(i));
    constante2=inv(gamam(i))^(-0.5*(H));
    densidad=constante1+constante2*exp((y-Xwn(:,1:1+d+m(i))*betam{i,1})'*(y-Xwn(:,1:1+d+m(i))*betam{i,1}));
    SUMA=SUMA+densidad;
    clear densidad;
end
MC=SUMA/M;
%
%   --  FIN --