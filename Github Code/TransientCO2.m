function [U,t,PL,xstar,P,M] = TransientCO2
% load parameters
dat = load('ModelParameterV1.0.mat'); 
DAY=load('AVGday.mat'); 
                                    
% declare the other dataset
ParLiana=dat.ParLiana;
ParTree=dat.ParTree;
Dem=dat.Dem;
DAY=DAY.DAY;

% declare variables
N = 3000;
x = linspace(Dem.x0,500,N).';   % vector with tree sizes
dx = x(2)-x(1);             % sizestep 
flip = N:-1:1;              % to flip the vector

% standard growth rates for liana and trees at current conditons
gL0=Dem.gL;    gT0=Dem.gc;  
trans = Dem.xstar/Dem.xstarL;       % tranformation trick to set lianas and trees in same model dimension

% run the model to obtain the standard fecundity rates
[~, ~,FT0, FL0] = Demography(Dem, ParTree, ParLiana,'indirect', x, 'mod1');

%% --> DEFINE THE CO2 SCENARIO
    % CO2 increase scenario
    fit=load('CO2_transient.mat');  % fitted curve to observed and predicted 
    CAincrease=fit.SSP2_fit(1975:2024).';  % CO2 increase between defined periods
    CAdiff=fit.SSP2_fit(2024)./fit.SSP2_fit(1975);  % difference in CO2 between start and end scenario
    % define the timeframe of the CO2 scenario
    tspin = 1000;    tchange=length(CAincrease);  trecover=500;
    tspan = 0:10:(tspin+tchange+trecover);          

%% --> PHYSIOLOGICAL MODEL TO DEFINE G UNDER ALTERED CO2 CONDITIONS
    % recalculate the corresponding carbon acquisition rates
    G1.L=Growth(ParLiana, DAY.DSL, DAY, CAincrease(1));  % Liana growth rate at start of scenario 
    G1.T=Growth(ParTree, DAY.DSL, DAY, CAincrease(1));   % Tree growth rate at start of scenario 
    G2.L=Growth(ParLiana, DAY.DSL, DAY, CAincrease(end));% Liana growth rate at the end of the scenario   
    G2.T=Growth(ParTree, DAY.DSL, DAY, CAincrease(end)); % Tree growth rate at the end of the scenario   
    G0.L=Growth(ParLiana, DAY.DSL, DAY, 400);            % Standard liana growth rate at current conditions    
    G0.T=Growth(ParTree, DAY.DSL, DAY, 400);             % Standard tree growth rate at current conditions 

%% --> BENCHMARK FOR THE STARTING CONDITIONS (i.e., in 1975) 
    % callibrated values at the start of the simulation
    betaL=0.41;     betaT=1;                             % beta of lianas and trees (the latter is standard)
    Dem.Dstar=25;   Dem.Xstar=D2x(Dem.Dstar,ParTree);    % forest height                
    Dem.PL=0.17;                                         % initial liana prevalence

%% --> TRANSLATE NEW G INTO CORRESPONDING GROWTH RATES
    % translate G to growth itself
    hetaT = (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
            (ParTree.phi).*(G2.T-G1.T)./Dem.gc./(CAdiff.^betaT-1);
    hetaL = (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
            (ParLiana.phi).*(G2.L-G1.L)./Dem.gL./(CAdiff.^betaL-1);
    % save the standard growth rates as to compare to
    gL1= gL0 + (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
        (ParLiana.phi*hetaL).*(G1.L-G0.L);  % Liana growth at start of CO2 scenario
    gc1= gT0 + (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
        (ParTree.phi*hetaT).*(G1.T-G0.T);  % Tree growth at start of CO2 scenario
    gL2= gL0 + (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
        (ParLiana.phi*hetaL).*(G2.L-G0.L);  % Liana growth at end of CO2 scenario
    gc2= gT0 + (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./...
        (ParTree.phi*hetaT).*(G2.T-G0.T);  % Tree growth at end of CO2 scenario
    Dem.gc=gc1;    Dem.gL=gL1;      % save the growth rates at the start of the scenario (i.e., 1975)


%% --> RUN THE MODEL WITH CHANGING GROWTH AND PLOT

    % run model for the start conditions (i.e., 1975) from where our simulation
    % will start from and save the obtained fecundity rates
    [PL1, ~,FT1, FL1, Dens] = Demography(Dem, ParTree, ParLiana,'indirect', x, 'mod1');
    Dem.FT=FT1; Dem.FL=FL1;

    % calculate the number of individual trees and understory trees at the
    % initial state of the run (i.e., 1975)
    guL1 = gL1*1.2*trans;  
    U1(1:N,1) = Dens.Nu+Dens.Nc;
    U1(N+1:2*N,1)=Dens.Ni;
    U1(2*N+1:3*N,1)=FL1*PL1/guL1*exp(-Dem.muL./guL1*x); 
      
    % run the model step by step, where initial conditions are provided as
    % vector U1 (as define just above). The model provides per simulated
    % timestep the forest structure in infested, non infested, and
    % liana recruits.
    b=nan(2,1);
    b(1)=(gc2./gc1-1)./length(CAincrease);  
    b(2)=(gL2./gL1-1)./length(CAincrease);
    opts = odeset('RelTol',1e-6,'AbsTol',1e-6);     %,'Jacobian',@Fjac);

    [t,U] = ode78(@(t,y) fun0(t,y), tspan, U1, opts);


%% PLOTTING OF THE FIGURE
    % define the CO2 scenario profile to plot
    CA1=ones(1,tspin)*CAincrease(1);
    CA2=CAincrease;
    CA3=ones(1,trecover)*CAincrease(end);
    CAperiod=[CA1(:); CA2(:); CA3(:)];

    % here we make a  selection of the timeframe plotted in the manuscript
    i1=find(tspan==900);  i2=find(tspan==1500);

    % Define empty matrixes for liana prevalence and tree structure
    AL = zeros(size(tspan));
    xstar = zeros(size(tspan));

    % for the requested timeframe we want to plot, we calculate the
    % distribution of infested, non infested, and recruiting lianas
    for j=i1:i2
        % extract the distribution of canopy trees and liana recruits
        LF(:,1) = U(j,1:N);           % liana free trees
        LI(:,1) = U(j,N+1:2*N);       % liana infested trees
        LL(:,1) = U(j,2*N+1:3*N);     % liana recruits
        % calculate canropy area
        Atree = -cumtrapz(x(flip),(LF(flip)+LI(flip)).*x(flip).^(1/(Dem.theta-1)));
        und = Atree(flip)>1;          % allocate who is undersatory
        can = ~und;                   % allocate who is canopy
        AL(j) = trapz(x(can),LI(can).*x(can).^(1/(Dem.theta-1)));   % calcualte the area occupied by lianas
        xstar(j)=x(find(can,1));      % finds the smallest diameter size belonging to the canopy     
    end

%% --> PLOT FIGURES
  figure(3)

    % plot the CO2 scenario
    subplot(311)
    area([tspin tspin+tchange], [450 450]); hold on;     % plot area in which CO2 increases
    plot(CAperiod,'k','LineWidth',2); 
        xlim([900 1500]);
        ylim([320 440]);
        ylabel('CO_2 (ppm)');
        xticks([1025 1125 1225 1325 1425]);
        xticklabels({'2000','2100','2200', '2300', '2400'});
        xlabel('year (C.E.)');
    
    % plot area infested by lianas  
    subplot(312)
    area([tspin tspin+tchange], [450 450]); hold on;     % plot area in which CO2 increases
    plot(tspan,AL*100,'k', 'LineWidth',2);
        xlim([900 1500]);
        ylim([17 23]);
        ylabel('P_L (%)');
        xticks([1025 1125 1225 1325 1425]);
        xticklabels({'2000','2100','2200', '2300', '2400'});
        xlabel('year (C.E.)');

     % plot the forest structure (transformed to the D dimension with x2D
     % function)
     subplot(313)
     area([tspin tspin+tchange], [450 450]); hold on;     % plot area in which CO2 increases
     plot(tspan,x2D(xstar, ParTree),'k', 'LineWidth',2);
        xlim([900 1500]);
        ylim([20 50]);
        ylabel('D* (cm)');
        xticks([1025 1125 1225 1325 1425]);
        xticklabels({'2000','2100','2200', '2300', '2400'})
        xlabel('year (C.E.)');




%% DEFINING THE TRANSIENT MODEL FUNCTIONS
  function dU = fun0(t,U)
      % forest demography
      U(U<0)=0;U=real(U);
      LF = U(1:N,1);        % liana free trees
      LI = U(N+1:2*N,1);    % liana infested trees
      LL = U(2*N+1:3*N,1);  % liana recruits
      % calculate canopy area
      Atree = -cumtrapz(x(flip),(LF(flip)+LI(flip)).*x(flip).^(1/(Dem.theta-1)));
      und = Atree(flip)>1; % which trees are undersatory
      can = ~und;          % these trees are canopy
      i = find(can,1);     % find the index of the smallest tree belonging to the caonpy 
      [t,i]                % to return how our progress of the run is going

      % define the corresponding growth rates based on the timing of the 
      % run (i.e., linked to the CO2 experienced at that moment)
        if t<=tspin
            gc = gc1;
            gL = gL1;
        elseif t>tspin && t<tspin+tchange
            gc = gc1*(1 + b(1)*(t-tspin));
            gL = gL1*(1 + b(2)*(t-tspin));
        elseif t >= tspin+tchange
            gc = gc2;
            gL = gL2;
        end
        
        % adjust the fecundity rates accordingly
        FL = FL1 *gL/gL1;
        FT = FT1 *gc/gc1;

        % define the distict growth rates for the different cathegories of
        % trees
        gu = gc/2;                                 % understory growth rate set at half canopy growth rate
        guL = gL*1.2*trans;                        % understory liana growth rate, but transformed to match that of trees 
                                                   % (trick to ease the modelling definitions)
        gLF = (gc-gu)./(1+exp(-(x-x(i)))) + gu;    % for Trees Free
        gLI = (0 -gu)./(1+exp(-(x-x(i)))) + gu;    % for Trees Infested

        % calculate the liana prevalence and defining the distinct
        % mortality rates
        PL = trapz(x(can),LI(can).*x(can).^(1/(Dem.theta-1)));
        mLF = ones(N,1)*Dem.mu;  mLF(can)=Dem.mc;
        mLI = ones(N,1)*Dem.mu;  mLI(can)=Dem.mi;

    % CALCULATION OF THE LIANA INFESTATION RATE
    % calculate lambda and use the RL obtained from the just declared liana
    % numbers arrived at xstar
        lambda = zeros(N,1);    
        x0L=x(i)./trans;       % size of liana to reach canopy
        lambda(can) = gL./ log( (x(can)./x0L).^(1/(Dem.theta-1))./guL./LL(i)./gL -1);     

    % NEW FOREST CONDITIONS
    dU = nan(3*N,1);
        % for trees
        u1 = zeros(N+1,1);
        u1(1) = FT*(1-PL)/gu;
        u1(2:N+1) = LF;
        dgdx = (gLF-gu).*(gc-gLF)/(gc-gu);
        dU(1:N,1) = -(mLF+dgdx+lambda).*LF - gLF.*diff(u1)/dx + Dem.xi*LI;

        u1(1) = 0;      % set initial value to 0
        u1(2:N+1) = LI; 
        dgdx = (gLI-gu).*gLI/gu;
        dU(N+1:2*N,1) = lambda.*LF - (mLI+dgdx+Dem.xi).*LI - gLI.*diff(u1)/dx;

        % calculate the number of lianas per size cohort
        u2 = zeros(N+1,1);
        u2(1) = FL*PL/guL;
        u2(2:N+1) = LL;
        dU(2*N+1:3*N,1) = -Dem.muL.*LL - guL.*diff(u2)/dx;

  end
end
