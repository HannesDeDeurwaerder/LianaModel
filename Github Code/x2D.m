function D = x2D(x,par)

%% compute diamater from x plant dimention

D = (x.^(1/(par.theta-1))/par.beta).^(1/par.c);
% x = (par.beta*D.^par.c).^(par.theta-1);