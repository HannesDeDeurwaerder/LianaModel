
function G = Growth(PAR, Q, DAY, Ca)

cf=2*12e-9*3600*24*365;      % conversion factor from umol CO2 m-2 s-1 to Kg dry mass yr-1

for i = 1:length(DAY.Time)
% DRY season
    dat.D= struct('I',DAY.D.I0(i),'Ca',Ca,'TL',DAY.D.Tc(i),'D',DAY.D.D(i),'psiS',DAY.D.PsiS(i),'H',25,'LightLimited','y');
    out.D=An_StomatOpt(PAR.Vcmax25, PAR.Kmax, PAR.p50, PAR.gamma, dat.D, PAR.cost);
            Gstar.D(i)=out.D.LL.Gstar;
             
% WET season      
    dat.W= struct('I',DAY.W.I0(i),'Ca',Ca,'TL',DAY.W.Tc(i),'D',DAY.W.D(i),'psiS',DAY.W.PsiS(i),'H',25,'LightLimited','y');
    out.W=An_StomatOpt(PAR.Vcmax25, PAR.Kmax, PAR.p50, PAR.gamma, dat.W, PAR.cost);
            Gstar.W(i)=out.W.LL.Gstar;                
end

G=(Q.*mean(Gstar.D) + (1-Q).*mean(Gstar.W))*cf; 
end