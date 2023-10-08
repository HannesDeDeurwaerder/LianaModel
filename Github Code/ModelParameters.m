% Parameter declaration
%--------------------------------------------------------------------------
% Author: Hannes De Deurwaerder, Matteo Detto, Stefan Schnitzer, Marco Visser, Stephen Pacala
% created: 2023/04
% Copyright: (c) H.P.T. De Deurwaerder et al.
% Email: Hannes_de_deurwaerder@hotmail.com
% 
% Description:
% This script will declare all used parameters


%% model parameters - overall params
    % graphical parameters
    col.L      = [0.9882    0.5804    0.1098];          % graphical parameters
    col.T      = [0.1451    0.6784    0.0627];          % graphical parameters

%% model parameters - Trees Parameter set

    % trees - physiology model
    ParTree.Pleaf=              2/3;               % how much more sensitive is the leaf compared to root and stem in P50 response (Scofoni and Sack et al.)
    ParTree.Vcmax25=            39.9431;           % [umol m-2 s-1]   from Norby et al 2017
    ParTree.KmaxStem =          5.9650;            % [in mmol m-2 s-1 MPa-1] 
    ParTree.Kmax(1,:) =         ParTree.KmaxStem;  % [in mol m-1 s-1 MPa-1]          % roots
    ParTree.Kmax(2,:) =         ParTree.KmaxStem;  % [in mol m-1 s-1 MPa-1]          % stem (via database)
    ParTree.Kmax(3,:) =         ParTree.KmaxStem;  % [in mol m-1 s-1 MPa-1]          % leaves
    ParTree.p50(1,:) =          0;                 % [MPa]                           % roots
    ParTree.p50(2,:) =          -1.7;              % [MPa]                           % stem (via database)
    ParTree.p50(3,:) =          -1.7*ParTree.Pleaf;% [MPa]                           % leaves
    ParTree.gamma(1,:) =        2;                 % the slope of the vulnerability curve 
    ParTree.gamma(2,:) =        2;                 % the slope of the vulnerability curve 
    ParTree.gamma(3,:) =        2;                 % the slope of the vulnerability curve 

    % trees cost function
    ParTree.t1 =                2;                 % exponent of cost function
    ParTree.TLP =               -2;                % Turgor loss point, obtained from Marechaux et al, 2017
    ParTree.c0 =                2.6072;            % cost on cavitation
    ParTree.cost =              @(x) ParTree.c0*(x/ParTree.TLP).^ParTree.t1;  % Costfunction

    % Trees allometry
    ParTree.beta =              0.66;              % factor of canopy allometry --> obtained from Martinez-Cano et al. 2019
    ParTree.psi =               0.0933;            % factor for biomass allometry --> obtained from Martinez-Cano et al. 2019
    ParTree.c =                 1.34;              % exponent of canopy allometry--> obtained from Martinez-Cano et al. 2019 
    ParTree.d =                 2.56;              % exponent of biomass function --> obtained from Martinez-Cano et al. 2019 
    ParTree.eta =               10;                 % metabolic costs related to nutrient acquisition, cost of biomass construction, cost of NSC storage,... (OPTIMIZED)
    ParTree.fi =                0.2;               % how much of G allocated to reproduction [between 0 -1]
    ParTree.theta =             ParTree.d/ParTree.c;                    % expoenent biomass~crown area
    ParTree.phi =               ParTree.psi/ParTree.beta^ParTree.theta; % factor biomass~crown area


%% model parameters - Lianas Parameter set

    % Lianas - physiology model
    ParLiana.Pleaf =            2/3;                % how much more sensitive is the leaf compared to root and stem in P50 response (Scofoni and Sack et al.)
    ParLiana.Vcmax25 =          50.0385;            % [umol m-2 s-1]   from Norby et al 2017 
    ParLiana.KmaxStem =         8.0828;             % [in micromol m-2 s-1 MPa-1] 
    ParLiana.Kmax(1,:) =        ParLiana.KmaxStem;  % [in mol m-1 s-1 MPa-1]         % roots
    ParLiana.Kmax(2,:) =        ParLiana.KmaxStem;  % [in mol m-1 s-1 MPa-1]         % stem (via database)
    ParLiana.Kmax(3,:) =        ParLiana.KmaxStem;  % [in mol m-1 s-1 MPa-1]         % leaves
    ParLiana.p50(1,:) =         0;                  % [MPa]                          % roots
    ParLiana.p50(2,:) =         -1.16;               % [MPa]                         % stem (via database)
    ParLiana.p50(3,:) =         -1.16*ParLiana.Pleaf;% [MPa]                         % leaves
    ParLiana.gamma(1,:) =       2;                  %  the slope of the vulnerability curve 
    ParLiana.gamma(2,:) =       2;                  %  the slope of the vulnerability curve 
    ParLiana.gamma(3,:) =       2;                  %  the slope of the vulnerability curve 

    % Lianas - cost function
    ParLiana.t1 =               2;                  % exponent of cost function
    ParLiana.TLP =              -2;                 % Turgor loss point -->  obtained from Marechaux et al, 2017 
    ParLiana.c0 =               4.7664;             % cost of cavitation --> orriginally this was 4
    ParLiana.cost =             @(x) ParLiana.c0*(x/ParLiana.TLP).^ParLiana.t1;  % Costfunction

    % Lianas allometry
    ParLiana.beta =             0.4123;             % factor of canopy allometry --> obtained from caclulation from BCI survey
    ParLiana.psi =              0.2190;             % factor for biomass allometry --> obtained from Schnitzer et al. 2006
    ParLiana.c =                1.4039;             % exponent of canopy allometry--> obtained assuming theta is the same for tree
    ParLiana.d =                2.682;              % exponent of biomass function --> obtained from Schnitzer et al. 2006
    ParLiana.eta =              10;                 % metabolic costs related to nutrient acquisition, cost of biomass construction, cost of NSC storage,... (OPTIMIZED)
    ParLiana.fi =               0.2;                % how much of G allocated to reproduction [between 0 -1]
    ParLiana.theta =            ParLiana.d/ParLiana.c;         % expoenent biomass~crown area
    ParLiana.phi =              ParLiana.psi/ParLiana.beta^ParLiana.theta;  % factor biomass~crown area

%% demographic model parameters
    Dem.PL =                    0.2;                    % fraction of liana occupancy
    Dem.fr =                    0.2;                    % fraction of liana occupancy
    Dem.mct =                   0.021;                  % total of all canopy trees (i.e. both infested as non infested trees)------>  obtained from BCI analysis
    Dem.mu =                    Dem.mct;                % mortality term understory tree ---> see BCI analysis
    Dem.v =                     0.021;                  % lethality from Visser et al. 2018
    Dem.mc =                    Dem.mct-(Dem.fr*0.021); %mortality of canopy trees --> (calculated given Dem.fr of canopy is liana infested and knowing that Dem.mi=Dem.mc+0.021)      
    Dem.mi =                    Dem.mc+Dem.v;           % mortality term liana infested tree --> obtained from Visser et al. (2018)
    Dem.muL =                   0.0709;                 % from Schintzer 2021 (smallest liana group)

    Dem.D0 =                    5;                      % minimal diameter considered in the model (in cm)
    Dem.x0 =                    D2x(Dem.D0,ParTree);    % D0 transformed in the x dimension
    Dem.Dstar =                 33;                     % obtained from the BCI analysis
    Dem.xstar =                 D2x(Dem.Dstar,ParTree); % Dstar transformed in the x dimension
    Dem.DstarL =                1.6;                    % minimal diameter for a liana to reach the canopy
    Dem.xstarL =                D2x(Dem.DstarL,ParLiana); %DstarL transformed in the x dimension

    Dem.theta =                 ParLiana.d/ParLiana.c;                     % expoenent biomass~crown area (same for liana and tree)
    Dem.phiL  =                 ParLiana.psi/ParLiana.beta^ParLiana.theta; % factor biomass~crown area liana
    Dem.phiT  =                 ParTree.psi/ParTree.beta^ParTree.theta;    % factor biomass~crown area tree

    Dem.xi =                    0.002;                   % Liana shedding, extracted from Visser et al, 2018
    Dem.gct =                   0.7143;                  % growth rate in ther x-dimension. --> Obtained from the BCI analysis (will eventually used to parameterize the physiological model)
    Dem.gc=                     Dem.gct/(1-Dem.fr);      % growth rate at the X-level for canopy species (hence accounting for higher mortality of total canopy due to liana infested trees included)
    Dem.gu =                    Dem.gc/2;                % growth rate at the X-level for the understory individuals, obtained via the BCI analysis where we saw it was half of the uper canopy
    Dem.gL =                    ParLiana.beta^(ParLiana.theta-1)*ParLiana.c*(ParLiana.theta-1)*...
                                Dem.DstarL^(ParLiana.c*(ParLiana.theta-1)-1)*0.05;    % growth rate at the X-level for lianas (growth in D level of liana = 0.05cm y-1 for 1 cm Liana according to the Liana census at BCI by Visser and Schnitzer)


  %% write away all parameter sets  
    save('ModelParameterV1.0.mat','col','ParLiana','ParTree','Dem')

