function output = An_StomatOpt(Vcmax25,k,p50,as, dat, CostFun)
% THIS FUNCTION IS OUR PHYSIOLOGICAL MODULE AS DESCRIBED IN THE MAIN MANUSCRIPT
% Function: Net Photosyntesis with stomata optimization coupled with a 3-element hydraulic system

%% INPUTS
%   Vcmax25:    maximum carboxillation velocity at reference temperature [mumol/m2/s]
%   k:          vector [kr kx kL] max hydraulic conductance per unit of leaf area
%               for root, xylem and leaf [mmol m-2 s-1 Mpa-1]
%   p50:        vector [0 p50x p50L]  water potential at 50% loss of conductivity
%               for xylem and leaf (not defined for root) [Mpa]
%   as:         vector [a asx asl] inflection of vulnerability curve for
%               xylem and leaf (not defined root) [Mpa-1]
%   c0:         cost parameter associated with psiL [mumol CO2 / Mpa-t / m2]
%   t:          nonlinear cost exponent of cost function
%
% dat: meteorological inputs structure
%       dat.TL:     leaf temperature (C)
%       dat.D:      VPD  (Kpa)
%       dat.I:      incident PAR (mumol m-2) (it can be a vector)
%       dat.Ca:     ambient CO2 concentration (ppm)
%       dat.H:      height of the leaf (m)

%% OUTPUTS
%   gs:             Stomatal conductance [mmol m-2 s-1]
%   An:             Net photosynthetic rate [mumol CO2 m-2 s-1]
%   E:              Evapotranspiration [...]
%   Rd:             Dark leaf respiration [mumol CO2 m-2 s-1]
%   psi0:           water potential at stem base [MPa]
%   psix:           xylem water potential [MPa]
%   psiL:           Leaf water potential [MPa]
%   iWUE:           instrinsic water use efficiency [...]
%   meta:           information on limitation status


%% EXAMPLE:
% dat = struct('I',1000,'Ca',400,'TL',25,'D',1.5,'psiS',-0.5,'H',30);
% Vcmax25 = 40;
% Kmax = [5 10 10];
% p50 = [0 -1.5 -1];
% as = [0 6 3];
% c0 = 1;
% t = 2;
% [An,Ac,Aj,gs0,Rd,psiL,psix,psi0]  = An_StomatOpt(Vcmax25,Kmax,p50,as,c0,t,dat);

%% ALGORITHM:
% The optimization is based on maximization of An - cost
% where water related cost = c0 * psiL^t
% psiL is computed from a three-element hydrualic model root-xylem-leaves
% soil water potential is the weighted average by root distribution
% For example, using Jackson 1996, psiS can be computed as
%       beta = 0.961; f=zeros(4,1);
%       f(1) = 1-exp(log(beta)*z1);
%       f(2) = exp(log(beta)*z2)-exp(log(beta)*z3);
%       f(3) = exp(log(beta)*z3)-exp(log(beta)*z4);
%       f(4) = exp(log(beta)*z4);
%       psiS = psi*f
% where psi is a vector of water potential at differnt depths z1, z2, z3 and z4 
% Photosyntesis is computed from the Farqhuar model with temperature
% dependent kynetics from Medelyn (2002)

%% PARAMETERIZATION
Tg =            25;                 % plant growing temperature
R =             8.414;              % J mol-1 K-1  
Kc25 =          404.9;              % mubar
E.Kc =          79430;              % J mol-1   
Ko25 =          278.4;              % mbar
E.Ko =          36380;              % J mol-1 
Gstar25 =       42.75;              % mumol/mol
E.Gstar =       37830;              % mubar
Rd25 =          0.015*Vcmax25;      % mumol m-2 s-1
E.Rd =          46390;              % J mol-1 

% Vcmax (from Medelyn et al., 2002)
Ha.Vcmax =      66560;              % J mol-1  
Hd.Vcmax =      200000;             % J mol-1  
S.Vcmax =       638-1.07*Tg;        % J mol-1 K-1 Kattge and Knorr (2007), table 3

% Jmax (from Medelyn et al., 2002)
Jmax25 =        1.48*Vcmax25+6.30;  % Norby et al (2017) for Panama
Ha.Jmax =       43965;              % J mol-1  
Hd.Jmax =       200000;             % J mol-1  
S.Jmax =        659-0.75*Tg;        % J mol-1  Kattge and Knorr (2007), table 3
theta =         .7;                 % unitless
alpha =         0.36;               % unitless; realized quantum yield

% environmental conditions
Ca =            dat.Ca;             % mubar
O  =            210;                % mbar
D  =            dat.D*10;           % mbar
TK =            dat.TL+273;         % Kelvin
psiS =          dat.psiS;           % Mpa
Hi =            dat.H;              % m

% temperature functions
Kc  =       Kc25*exp((TK-298)./(R*TK*298)*E.Kc);
Ko  =       Ko25*exp((TK-298)./(R*TK*298)*E.Ko);
Gstar  =    Gstar25*exp((TK-298)./(R*TK*298)*E.Gstar);
Rd  =       Rd25*exp((TK-298)./(R*TK*298)*E.Rd);

Vcmax  =    Vcmax25*exp((TK-298)./(R*TK*298)*Ha.Vcmax)...
                    .*(1+exp((298*S.Vcmax-Hd.Vcmax)/(298*R)))...
                    ./(1+exp((S.Vcmax*TK-Hd.Vcmax)./(R*TK)));
   
Jmax   =    Jmax25*exp((TK-298)./(R*TK*298)*Ha.Jmax)...
                    .*(1+exp((298*S.Jmax-Hd.Jmax)/(298*R)))...
                    ./(1+exp((S.Jmax*TK-Hd.Jmax)./(R*TK)));
      
Km =        Kc.*(1+O./Ko); 


% options for optimizer
opts = optimoptions(@fmincon,'Algorithm','active-set','Display','off',...
    'SpecifyConstraintGradient',false);

% for hydraulics
rho =     1e3;               % density of water [kg m-3]                            
g =       9.81;              % gravity [m s-2]
pz =      rho*g*Hi*1e-6;     % gravitational term with correction term from Pa to MPa

%% CARBON LIMITATION
% Farquahar
a  =  (Vcmax-Rd)/2;
b  =  (Ca+Km)/2;
c  =  Rd./2.*(Ca+Km) + Vcmax./2.*(2*Gstar-Ca+Km);

gc = fmincon(@hydraulics,0.001,-1,0,[],[],0,[],@nonlcon,opts);

Ac(:,1) = a + b.*gc - sqrt(b.^2.*gc.^2+c.*gc+a.^2);

% Evapotranspiration
E = 1.6*gc*D;

% Hydraulics
psi0 = psiS - E/k(1);
psix = p50(2) - pz + 1/as(2)*log((exp(as(2)*(psi0-   p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
psiL = p50(3) - pz + 1/as(3)*log((exp(as(3)*(psix+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);
cost = CostFun(psiL);

% Intrinsic WUE
dAndgs=b-((2.*b^2.*gc +c)/(2.*(sqrt(b^2.*gc^2+c.*gc+a^2))));
iWUE=dAndgs./dat.D;

% write output file for carbon limited conditions
output.CL.gs = 1.6*gc;
output.CL.An = Ac; 
output.CL.E = E;
output.CL.psi0 = psi0;
output.CL.psix = psix;
output.CL.psiL = psiL;
output.CL.cost = cost;
output.CL.Rd = Rd;       
output.CL.iWUE = iWUE;
output.CL.meta='carbon limited';
output.CL.Gstar = Ac - cost;

%% LIGHT LIMITATION
if strcmpi(dat.LightLimited,'y')  

    % Farquahar   
    J  =  (alpha.*dat.I + Jmax - sqrt((alpha.*dat.I+Jmax).^2 - 4*alpha*Jmax.*theta.*dat.I))./(2*theta);
    a  =  J/8-Rd/2;
    b  =  Ca/2+Gstar;
    c  =  Rd./2.*(Ca  + 2*Gstar) + J./2.*(Gstar - Ca/4);
    
    gj = fmincon(@hydraulics,0.001,-1,0,[],[],0,[],@nonlcon,opts);
    Aj = a + b.*gj - sqrt(b.^2.*gj.^2+c.*gj+a.^2);
    
    [An,index] = min([Ac Aj]);
    if index==1
        gs = gc;
    else
        gs = gj;
    end
    
    % Evaporation
    E = 1.6*gj*D;
    
    % hydraulics
    psi0 = psiS - E/k(1);
    psix = p50(2) - pz + 1/as(2)*log((exp(as(2)*(psi0-   p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
    psiL = p50(3) - pz + 1/as(3)*log((exp(as(3)*(psix+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);
    cost = CostFun(psiL);

    % Intrinsic WUE
    dAndgs=b-((2.*b^2.*gj +c)/(2.*(sqrt(b^2.*gj^2+c.*gj+a^2))));
    iWUE=dAndgs./dat.D;
    
    % write output file for light limited conditions
    output.LL.gs = 1.6*gj;
    output.LL.An = Aj; 
    output.LL.E = E;
    output.LL.psi0 = psi0;
    output.LL.psix = psix;
    output.LL.psiL = psiL;
    output.LL.cost = cost;
    output.LL.Rd = Rd;                                       
    output.LL.iWUE = iWUE;
    output.LL.meta='Light limited';
    output.LL.Gstar = Aj - cost;
else
    An = Ac;
    gs = gc;
end


%% MINIMIZING CONDITIONS (optimizer)
% Evapotranspiration
E = 1.6*gs*D;

%hydraulics
psi0 = psiS - E/k(1);
psix = p50(2) - pz + 1/as(2)*log((exp(as(2)*(psi0-   p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
psiL = p50(3) - pz + 1/as(3)*log((exp(as(3)*(psix+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);
cost = CostFun(psiL);

 % Intrinsic WUE
dAndgs=b-((2.*b^2.*gs +c)/(2.*(sqrt(b^2.*gs^2+c.*gs+a^2))));
iWUE=dAndgs./dat.D;
    
% write output file for mimimum conditions
output.M0.gs = 1.6*gs;
output.M0.An = An; 
output.M0.E = E;
output.M0.psi0 = psi0;
output.M0.psix = psix;
output.M0.psiL = psiL;
output.M0.cost = cost;
output.M0.Rd = Rd;                                         
output.M0.iWUE = iWUE;
output.M0.Gstar = An - cost;
output.M0.meta='minimum between Ac and Aj';


%% FUNCTIONS FOR SOLVER
function [y,p0,p1,p2] = hydraulics(x)

E = 1.6*x*D;

p0 = psiS - E/k(1);
p1 = p50(2) - pz + 1/as(2)*log((exp(as(2)*(p0   -p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);
p2 = p50(3) - pz + 1/as(3)*log((exp(as(3)*(p1+pz-p50(3))) + 1).*exp(-as(3)*E/k(3)) - 1);

y = CostFun(p2)  - a - b.*x + sqrt(b.^2.*x.^2+c.*x+a.^2);

end


function [c,ceq] = nonlcon(x)

E = 1.6*x*D;

p0 = psiS - E/k(1);
p1 = p50(2) - pz + 1/as(2)*log((exp(as(2)*(p0-p50(2))) + 1).*exp(-as(2)*E/k(2)) - 1);

c = 1 - (exp(as(3)*(p1+pz-p50(3))) + 1).*exp(-as(3)*E/k(3));
ceq = [];

end


end
