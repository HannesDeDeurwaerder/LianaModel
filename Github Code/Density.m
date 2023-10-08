function Dens = Density(param,dat, D,name)

RP = param(2);
y= Demography(param,dat);
q1=(y(1)+1)/ dat.FH;
q2=(y(2)+1)/((dat.FP/RP-1).*dat.FH);

% optain occupancy values
Ah= q1/(q1+q2)*100;
Ap = q2/(q1+q2)*100;
Dstar=param(1).^(1/(dat.d-dat.c));   %back to DBH (cm)
RP=param(2);

%-----------------------------------------------
% obtain the abundance per tree cohort
%----------------------------------------------

% declare
gc=dat.gc;
gu=dat.gc/2;      
gp=dat.gp;
d=dat.d;
xi=dat.xi;
v=dat.v;
mc=dat.mc;
mu=dat.mu;
thetaP=dat.thetaP;
thetaH=dat.thetaH;
phiP=dat.phiP;
phiH=dat.phiH;
FH=dat.FH;
FP=dat.FP;
x0=dat.x0;
c=dat.c;

% declare vectors
Nu=nan(size(D));
Nc=nan(size(D));
Ni=nan(size(D));

% now calculate all the tree cohort specific densities (#/m-2)
RH=FH.*(1-RP./FP);
x=D.^(d-c);
Xstar=Dstar.^(d-c);       
trans= (d-c)*D.^(d-c-1);  
use = x<Xstar;

    %understory   
    if sum(use)>0
        Nu(use)=RH./gu.*exp(-mu./gu.*(x(use)-x0));
    end     
    
    % canopy
        c1 =1/(thetaP+1); 
        c2=(1+thetaP-thetaH)/(thetaP+1);
        c3=phiP/phiH;
        mn = mc+v+xi;
        b2 =(mc+v)./c2./mn.*(RP.*c1.*c3.*gp.^thetaP).^c1;
        der=mc.*(x(~use)-Xstar)+b2.*(x(~use).^c2-Xstar.^c2);
        Nc(~use)=RH./gc.*exp((-mu./gu.*(Xstar-x0))-(der./gc));
        
    % infested
        lambda=(RP.*c1.*c3.*gp.^thetaP./(x(~use).^thetaH)).^c1;
        Ni(~use)=(lambda./mn).*Nc(~use);
        
    % express per m-2 and write away
    if strcmp(name,'und')
      Dens=Nu.*trans./phiH;           % # m-2
    elseif strcmp(name,'can')
      Dens=Nc.*trans./phiH;           % # m-2
    elseif strcmp(name,'inf')
      Dens=Ni.*trans./phiH;           % # m-2
    end
      
end