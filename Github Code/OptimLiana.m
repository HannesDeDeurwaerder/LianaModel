% function for optimizer
function [w,f] = OptimLiana(y,params)

%declare variables
    PAR=params.PAR;
    datw=params.datw;
    datd=params.datd;
    BMwet=params.BMwet;
    BMdry=params.BMdry;
    Psi_pd=params.Psi_pd;

% change in dataset used in the physiological model
    PAR.Kmax(1:3,:)=y(1);
    PAR.Vcmax25=y(2);
    PAR.cost= @(x) y(3)*(x/PAR.TLP).^PAR.t1;
     
% Run the physiological model for the studied growthfrom during wet season
     out=An_StomatOpt(PAR.Vcmax25, PAR.Kmax, PAR.p50, PAR.gamma, datw, PAR.cost);
            f(1)=abs(out.CL.An-BMwet(1))/BMwet(1);                  %optimize for AN      
            f(2)=abs(out.CL.gs-BMwet(2))/BMwet(2);                  % optimize for gs
            f(3)=abs((out.CL.psiL-Psi_pd(1))-BMwet(3))/(-BMwet(3)); %optimize for psi diff
            
% Run the physiological model for the studied growthfrom during dry season
     out=An_StomatOpt(PAR.Vcmax25, PAR.Kmax, PAR.p50, PAR.gamma, datd, PAR.cost);
            f(4)=abs(out.CL.An-BMdry(1))/BMdry(1);                  %optimize for AN      
            f(5)=abs(out.CL.gs-BMdry(2))/BMdry(2);                  % optimize for gs
            f(6)=abs((out.CL.psiL-Psi_pd(2))-BMdry(3))/(-BMdry(3)); %optimize for psi diff

 % Provide weights to the different parameters we benchmark to           
     weigths=[1 5 1 1 3 1].'; 
     w = f*weigths;

end