% Modelo semiparametrico de regresion via bases radiales y onduletas
%
% Funciones:
%   
%   bwn     	-       funcion principal
%   bwn_demo    -       demostracion de 'bwn'
%   bwn_pred    -       calcula una muestra de la distribucion predictiva un paso adelante usando el modelo de onduletas radiales
%   bwn_verint  -       calculael estimador por importancia de la verosimilitud integrada del modelo
%   bwn_verint0 -       calcula el estimador de Monte Carlo de la verosilitud integrada del modelo de bases radiales de ondueltas.
%   bwndil      -       dilata aleatoriamente una columna de Xwn  con base en X (original)
%   bwnmat      -       genera la matriz de regresion de bases de onduletas con la funcion
%   bwnmv       -       calcula la pseudo log-maxima verosimilitud del modelo necesaria para el Factor de Bayes
%   bwnnac      -       añade una columna de base de onduletas indice 'k_prop'
%   bwntras     -       traslada aleatoriamenet una columna de Xwn  con base en X (original)
%   muest_sr    -       genera una muestra aleatoria sin reemplazo de tamaño k de un vector de n dimensiones
