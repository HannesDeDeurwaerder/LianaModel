function [PL,Dstar,FT, FL, Dens] = Demography(dat,ParTree,ParLiana,type, x, mod)


% declared globally
xstar0T = dat.xstar;
xstar0L = dat.xstarL;

gc = dat.gc;
gu = dat.gc/2;
gL = dat.gL;
guL = dat.gL*1.2*xstar0T/xstar0L;  % corresponds to 25 year time xstar/guL = time to reach canopy

xi = dat.xi;
mc = dat.mc;
mu = dat.mu;
mi = dat.mc+dat.v;
muL = dat.muL;

theta = dat.theta;
theta1 = theta-1;

x0 = dat.x0;

opts = optimoptions('lsqnonlin','display','off');


% this part estimates PL or Dstar -->  DIRECT method
if strcmp(type,'direct')

    init=[dat.xstar dat.PL];   
    lower=[20 1e-3];     
    upper=[90 1]; 

    % additional variables known to us
    FT = dat.FT;
    FL = dat.FL;

%%

    % OPTIMIZE DIRECT MODEL
    %-----------------------
    par = lsqnonlin(@(x) modelDirect(x),init,lower,upper, opts); %run model
    xstar = par(1);
    PL = par(2);
    Dstar = x2D(par(1),ParTree);


    % OPTIMIZE DIRECT MODEL
    %-----------------------
    elseif strcmp(type,'indirect')
    
    % additional variables known to us
    PL = dat.PL;
    Dstar = dat.Dstar;
    xstar=D2x(Dstar, ParTree);

     FLmax = 1e+6;
    
    
    % optimize indirect model (with separate stepS for FL and FT) 
    FL = lsqnonlin(@(x) modelIndirect1(x),100,0,FLmax, opts); %run model
    I1 = integral(@(x) fun1(x,xstar,FL.*PL),xstar,1e+6);
    FT = 1./I1;

end



% CHECK IF USER WANTS COI AND DENSITY
%------------------------------------
if nargout>4
    
    if ~exist('x','var')
        xmax=1000;
        xmin=dat.x0;
        x=linspace(xmin, xmax ,1000);    % in cm
    end

    dx=x(2)-x(1);   
    Dens.Nu=zeros(size(x));
    Dens.Nc=zeros(size(x));
    Dens.Ni=zeros(size(x));

    % now calculate all the tree cohort specific densities (#/m-2)
        und = x<xstar;
        can=~und;
        s=find(can,1);

        RT=FT.*(1-PL);
        RL=FL.*PL;
        g1 = @(x) (gc-gu)./(1+exp(-(x-xstar))) + gu;

        %understory   
        I1 = log(g1(x)/g1(x0));
        I2 = mu*((x-x0)/gc- ...
            (1/gu-1/gc)*...
            log((gc+gu*exp((xstar-x)))/(gc+gu*exp((xstar-x0)))));
        Dens.Nu(und) = RT./g1(x0)*exp(-(I1(und)+I2(und)));

        I1 = log(g1(xstar)/g1(x0));
        I2 = mu*((xstar-x0)/gc- ...
            (1/gu-1/gc)*...
            log((gc+gu)/(gc+gu*exp((xstar-x0)))));
        Nustar= RT./g1(x0)*exp(-(I1+I2));

        % canopy       
        dgdx = @(x) (g1(x)-gu).*(gc-g1(x))/(gc-gu);
        I3 = zeros(size(x));
        for j=s:length(x)
            I3(j) = integral(@(xp) (mc + dgdx(xp) + mi/(mi+xi)*lambda(xp, xstar, RL, mod))./g1(xp),xstar,x(j));
        end

        Dens.Nc(can)=Nustar*exp(-I3(can));
        Dens.Ni(can)=(lambda(x(can),xstar,RL, mod)./(mi+xi)).*Dens.Nc(can);
        

    % prevalence
        Dens.COI=Dens.Ni./(Dens.Nc+Dens.Ni)*100;
    
    % express per m-2 and write away
        Dens.D = x2D(x,ParTree);
        Dens.x = x;
        B1= integral(@(x) RT/g1(x0).*exp(-mu./gu.*(x-x0)).*x.^(theta/(theta-1)),x0,xstar);
        B2= integral(@(x) RT.*(fun3(x,xstar,RL)+fun4(x,xstar,RL)),xstar,1e6);

        Dens.B = ParTree.phi*(B1+B2)*1e+4*1e-3;
        Dens.P = trapz(x, (Dens.Nu+Dens.Nc).*(dat.phiT*theta/theta1.*x.^(1/theta1).*g1(x)))*10;
end

%% THE DIFFERENT MODELS

% DIRECT MODEL
function y = modelDirect(param)
        
        % paremeters
        xstar = param(1);  % xstar
        RL = param(2)*FL;  % recruitment rate

        % integral calculations
        I1 = integral(@(x) fun1(x,xstar,RL),xstar,1e+6);
        I2 = integral(@(x) fun2(x,xstar,RL),xstar,1e+6);
        % output
    
        y(1) = FT.*I1 - 1;
        y(2) = FT*(FL-RL)*I2 - RL;

end

function y = modelIndirect1(param)

        % paremeters 
        FL = param(1);
        RL = FL.*PL;

        % integral calculations
        I1 = integral(@(x) fun1(x,xstar,RL),xstar,1e+6);
        I2 = integral(@(x) fun2(x,xstar,RL),xstar,1e+6);
        
        % output
        y = I2./I1 - PL/(1-PL);

end




% FIRST INTEGRAL
function y = fun1(x,xstar,RL)
        

        g1 = @(x) (gc-gu)./(1+exp(-(x-xstar))) + gu;
        dgdx = @(x) (g1(x)-gu).*(gc-g1(x))/(gc-gu);
        
        K0 = integral(@(xp) (mu+dgdx(xp))./g1(xp),x0,xstar);
        K1 = zeros(size(x));        

        for i=1:length(x)
            K1(i) = integral(@(xp) (mc+dgdx(xp)+(mc+mi./(mi+xi))*lambda(xp, xstar, RL, mod))./g1(xp),xstar,x(i));
        end
        
        y =  1/g1(x0).*exp(-K0).*exp(-K1).*x.^(1/theta1);


end

% SECOND INTEGRAL
function y = fun2(x,xstar,RL)
      
        g1 = @(x) (gc-gu)./(1+exp(-(x-xstar))) + gu;
        dgdx = @(x) (g1(x)-gu).*(gc-g1(x))/(gc-gu);
        
        K0 = integral(@(xp) (mu+dgdx(xp))./g1(xp),x0,xstar);
        K1 = zeros(size(x));        

        for i=1:length(x)
            K1(i) = integral(@(xp) (mc+dgdx(xp)+(mc+mi./(mi+xi))*lambda(xp, xstar, RL, mod))./g1(xp),xstar,x(i));
        end
        
        y =  1/g1(x0).*exp(-K0).*lambda(x, xstar, RL, mod)./(mi+xi).*exp(-K1).*x.^(1/theta1);

end



% LAMBDA CALCULATION
function y = lambda(x, xstar, RL, mod)

if      strcmp(mod,'mod1')       % standard model (see manuscript, table 1)
        pc = exp(-muL/guL.*xstar);
        x0L = xstar*xstar0L/xstar0T;
        y = gL./ log((x./x0L).^(1/theta1)./(gL*RL.*pc) -1);

elseif  strcmp(mod,'mod2')       % Model 2 (see manuscript, table 1)
        pc = exp(-muL/guL.*x);
        x0L = 0;
        y = gL./ ((theta/theta1*gL*x.^(1/theta1)./RL./pc+ x0L.^(theta/theta1)).^(theta1/theta) - x0L);

elseif strcmp(mod,'mod3')        % Model 3(see manuscript, table 1)
        pc = exp(-muL/guL.*xstar);
        x0L = x*xstar0L/xstar0T;
        y = gL./ log((x./x0L).^(1/theta1)./(RL.*pc));

end


end

%% for biomass
% THIRD INTEGRAL
    function y = fun3(x,xstar,RL)

        g1 = @(x) (gc-gu)./(1+exp(-(x-xstar))) + gu;
        dgdx = @(x) (g1(x)-gu).*(gc-g1(x))/(gc-gu);
        
        K0 = integral(@(xp) (mu+dgdx(xp))./g1(xp),x0,xstar);
        K1 = zeros(size(x));        

        for i=1:length(x)
            K1(i) = integral(@(xp) (mc+dgdx(xp)+(mc+mi./(mi+xi))*lambda(xp, xstar, RL, mod))./g1(xp),xstar,x(i));
        end
        
        y =  1/g1(x0).*exp(-K0).*exp(-K1).*x.^(theta/theta1);

        

end

% FOURTH INTEGRAL
    function y = fun4(x,xstar,RL)
      
        g1 = @(x) (gc-gu)./(1+exp(-(x-xstar))) + gu;
        dgdx = @(x) (g1(x)-gu).*(gc-g1(x))/(gc-gu);
        
        K0 = integral(@(xp) (mu+dgdx(xp))./g1(xp),x0,xstar);
        K1 = zeros(size(x));        

        for i=1:length(x)
            K1(i) = integral(@(xp) (mc+dgdx(xp)+(mc+mi./(mi+xi))*lambda(xp, xstar, RL, mod))./g1(xp),xstar,x(i));
        end
        
        y =  1/g1(x0).*exp(-K0).*lambda(x, xstar, RL, mod)./(mi+xi).*exp(-K1).*x.^(theta/theta1);

end


end