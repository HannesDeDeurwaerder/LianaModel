function [x, dx_dD] = D2x(D,par)

%% compute x from diamater plant dimention

% D = (x.^(1/(par.theta-1))/par.beta).^(1/par.c);
x = (par.beta*D.^par.c).^(par.theta-1);

dx_dD=(par.theta-1)*par.c*par.beta.^(par.theta-1).*D.^(par.c.*(par.theta-1)-1);