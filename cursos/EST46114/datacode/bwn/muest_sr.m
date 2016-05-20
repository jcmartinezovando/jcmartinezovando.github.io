function [mue] = muest_sr(x,k)
% Toma una muestra aleatoria sin reemplazo de tamaño k de un vector de n dimensiones
%
%   NOTA: Utilizado en 'bwn.m'
%                   k << n
%   ENTRADAS:
%           x       -       vector de dimension 'nx1'
%           k       -       tamaño de la muestra
%
%   SALIDAS:
%           mue     -       muestra aleatorio de dimension k

n=max(size(x));
if k > n
    error('El tamaño de muestra deseado es mayor a la dimensión de los datos.')
end
% incializamos 'mue'
mue=zeros(k,1);
for i=1:k
    indx=ceil(rand*n);
    mue(i)=x(indx);
    x(indx)=[];
    n=n-1;
end
%
%       --  FIN --