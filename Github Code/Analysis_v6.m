% Liana model - Main analysis
%--------------------------------------------------------------------------
% Author: Hannes De Deurwaerder, Matteo Detto, Stefan Schnitzer, Marco Visser, Stephen Pacala
% created: 2023/04
% Copyright: (c) H.P.T. De Deurwaerder et al.
% Email: Hannes_de_deurwaerder@hotmail.com
% 
% Description:
% The following code was used for analysis and figure making
%
% Required datasets:
% The 50-ha tree census and the eddy covariance data for Barro Colorado
% Island (Panama) can be accessed via the DRYAD repository 
% (https://doi.org/10.15146/5xcp-0d46, and, 
% https://doi.org/10.5061/dryad.3tx95x6j5, respectively). 
% For the availability of the liana census and canopy occupancy index 
% data, readers are referred to Schnitzer et al., (2021) and Visser 
% et al., (2018) respectively.
%
% NOTE: loading of the base datasets will contain:
% 'Dem'     - All parameters related to the demography model
% 'ParLiana'- All liana related parameters
% 'ParTree' - All tree related parameters
% 'DAY'     - All environmental parameters for an average dry and wet day at BCI
%             Obtained from the flux tower data as described in the manuscript

%% Figure 3 
% Simulated and observed average liana canopy occupancy index (COI) with host size

% Clean and load required datasets
clear;   clc; 
run ModelParameters;           % source code with all model parameters
load('ModelParameterV1.0.mat');% loading all required model parameters
load('AVGday.mat');            % dataset of required environmental variables for an average 
                               % wet and dry day at BCI

% Upload template figure
% --> due to author restrictions of the datasets, we provide a figure
% template containing both the COI data and the tree demography data for
% BCI. When proof of access to the raw data is provided, we can provide 
% the processing code for COI and BCI demography upon request
TempFig=openfig("Fig3_Template.fig");
hold on


% -------------  Pannel A 

    % declare variables and graphical parameters
    x = linspace(Dem.x0,500,3000).';    % vector of different sizes to explore
    subplot('Position',[0.1411    0.4333    0.7750    0.5466])
    mod=["mod1";"mod2";"mod3"];
    LTY=[ "-"; ":"; "--"];
    LWD=[3,1.5,1.5];
  
    % run the model
    for q=1:length(mod)
    [PL,Dstar,FT, FL,  Dens] = Demography(Dem, ParTree, ParLiana,'indirect', x, mod(q));
    plot(Dens.D, Dens.COI, LTY(q), 'LineWidth',LWD(q), 'color','k');
    bb= find(~isnan(Dens.COI) == 1, 1);
    plot([Dstar Dstar],[0 Dens.COI(bb)],LTY(q), 'LineWidth',LWD(q), 'color','k');
    end
    plot([Dem.D0 Dstar],[0 0], 'LineWidth',3, 'color','k');

    % Figure cleaning
    h = get(gca,'children');
    legend(h([9 1 5 3]),'Data', 'Model 1', 'Model 2', 'Model 3', 'location','southwest');
    legend boxoff 
    set(gca, 'FontSize',12,'xscale','log')


% -------------  Pannel C

    % declare variables and bins
    Nbin=7;     Dlim=110;       Dmax=1000;
    bin=linspace(35,Dlim,Nbin-1);
    bin(1)=Dem.Dstar;
    bin(Nbin)=Dmax;
    Area_inf=nan(6,1);   % bin model data - infested trees
    Area_non=nan(6,1);   % bin model data - non-infested trees
    Dbin= mean([bin(1:end-1); bin(2:end)]);
    fbin = nan(length(bin)-1,1);

    % run the model
    x=linspace(Dem.xstar, D2x(Dmax, ParTree),3000);
    [PL, Dstar, FT, FL,  Dens] = Demography(Dem, ParTree, ParLiana,'indirect', x, 'mod1');
    % bin the model output
    for i=1:6      
        use=Dens.D>=bin(i) & Dens.D<bin(i+1);
        Area_inf(i)=trapz(Dens.Ni(use).*ParTree.beta.*Dens.D(use).^ParTree.c).*100;
        Area_non(i)=trapz(Dens.Nc(use).*ParTree.beta.*Dens.D(use).^ParTree.c).*100;
    end


%%% plot figure
    subplot('Position',[0.1354    0.1853    0.7750    0.1577])
    Dbin2=Dbin;Dbin2(end)=Dbin(end-1)+20;
    
    y=Dbin2+((Dbin2(2)-Dbin2(1))./5);
    b=bar(y,[Area_non Area_inf], 0.4, 'stacked'); 
    b(1).FaceColor=col.T;
    b(2).FaceColor=col.L;


%% Figure 4
%Changes in liana and tree growth rates, liana canopy prevalence (PL) 
% and plot biomass following an increase in ambient CO2 (CO2 increase; left panels) 
% and in dry season length (Drought shift; right panels).
clear;   clc;

% Load Base and declare global variables
run ModelParameters;
load('ModelParameterV1.0.mat'); 
load('AVGday.mat')  ;
load('Co2.mat'); 

% declare variables related to tested scenario
mod="mod1";                             % liana infestation model used
LTY=[ "-"; ":"; "--"];                  % linetypes for plotting

N=6;                                    % step sized
Qrange=linspace(1,1.1,N)*DAY.DSL;       % explored dry season length
CArange=linspace(400,425,N);            % explored CO2 range
x = linspace(Dem.x0,500,3000).';        % set of diameter ranges

g0L=Dem.gL;         % starting growth rate liana at current conditions BCI
g0T=Dem.gc;         % starting canopy growth rate trees at current conditions BCI

%------>>> DEFINE STANDARD MODEL AND PARAMETERS FOR COMPARISON

figure(2);

    % run standard model
    [PL,Dstar,FT, FL] = Demography(Dem, ParTree, ParLiana,'indirect', x, mod);
    Dem.FT=FT;Dem.FL=FL;            % save standard fecundity rates
    F0L=Dem.FL;   F0T=Dem.FT;       % save them as the standard starting values

    % altering parameters
    beta=[1.0 1.6 0.4];         % explored values for beta (relative growth response to increasing CO2 increase)
    hetaT = (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./(ParTree.phi)...
            .*(GTc(end)-GTc(1))./Dem.gc./((425/400).^beta-1);
    hetaL = (Dem.theta-1)./Dem.theta*(1-ParTree.fi)./(ParLiana.phi)...
            .*(GLc(end)-GLc(1))./Dem.gL./((425/400).^[beta(1) beta(3) beta(2)]-1);

%------>>> CO2 FERTILIZATION

    % run model for the various betas and simultatiously plot
    for p = 1:3
        % empty matrixes 
        PL=nan(N,1);
        Dstar=nan(N,1);
        B=nan(N,1);
        % calculate the new growth conditions in comparison with standard
        gL = g0L + (Dem.theta-1)/Dem.theta*(1-ParTree.fi)/(ParLiana.phi*hetaL(p))*(GLc-G0.L);
        gT = g0T + (Dem.theta-1)/Dem.theta*(1-ParTree.fi)/(ParTree.phi*hetaT(p))*(GTc-G0.T);
        % adjust the fecundity values accordingly
        FL = F0L *gL/g0L;
        FT = F0T *gT/g0T;

        parfor i = 1:length(gL)     
            [p i]
            % write away new values in new parameter set for the demographic model
            New=Dem;
            New.gL=gL(i);    New.gc=gT(i);     
            New.FL=FL(i);    New.FT=FT(i);
    
            % run the model with new conditions
            [PL(i),Dstar(i),~,~,Dens] = Demography(New, ParTree, ParLiana,'direct', x, mod);
            B(i)=Dens.B;        % extract the biomass for plotting
        end
    
    % make subplots
    subplot(421)
        % conversion from x to D dimension for plotting
        gD = gL/((Dem.theta-1).*1.6^(ParLiana.c*(Dem.theta-1)-1)...
                *ParLiana.beta^(Dem.theta-1)*ParLiana.c);   
            
        plot(CArange, gD,LTY(p),'color', col.L, 'LineWidth',2); hold on
        set(gca,'xticklabel',[]); 
        ylabel({'Diam. growth rate';'(cm yr^-^1)'}); 
        title('CO_2 increase');
        ylim([0.047 0.055]); xlim([CArange(1) CArange(end)]);
    
    subplot(423)
        % conversion from x to D dimension for plotting
        gD = gT/((Dem.theta-1).*33^(ParTree.c*(Dem.theta-1)-1)...
                 *ParTree.beta^(Dem.theta-1)*ParTree.c);
            
        plot(CArange, gD, LTY(p),'color', col.T, 'LineWidth',2); hold on
        set(gca,'xticklabel',[]);    
        ylabel({'Diam. growth rate';'(cm yr^-^1)'});     
        ylim([0.47 0.54]); xlim([CArange(1) CArange(end)]);
    
     subplot(425)
        plot(CArange, PL*100, LTY(p), 'LineWidth',2, 'color','k'); hold on;    
        set(gca,'xticklabel',[]);     
        ylabel({'Liana prevalence';' (P_L, %)'});
        ylim([17 27]); xlim([CArange(1) CArange(end)]);
    
    subplot(427)
        plot(CArange, B,LTY(p), 'LineWidth',2, 'color','k'); hold on
        xlabel('ambient CO_2 (ppm)');        
        ylabel({'Biomass';'(Mg ha^{-1})'})
        ylim([240 360]); xlim([CArange(1) CArange(end)]);
        pause(.1)
    end        

%------>>> DROUGHT INTENSIFICATION

    % run model for the various betas and simultatiously plot
    for p = 1:3
        % empty matrixes 
        PL=nan(N,1);
        Dstar=nan(N,1);
        B=nan(N,1);
        % calculate the new growth conditions in comparison with standard
        gL = g0L + (Dem.theta-1)/Dem.theta*(1-ParTree.fi)/(ParLiana.phi*hetaL(p))*(GL-G0.L);
        gT = g0T + (Dem.theta-1)/Dem.theta*(1-ParTree.fi)/(ParTree.phi*hetaT(p))*(GT-G0.T);
        % adjust the fecundity values accordingly
        FL = F0L *gL/g0L;
        FT = F0T *gT/g0T;

        parfor i = 1:length(gL)
            [p i]
            % write away new values in new parameter set for the demographic model
            New=Dem;
            New.gL=gL(i);    New.gc=gT(i);
            New.FL=FL(i);    New.FT=FT(i);

            % run the model with new conditions
            [PL(i),Dstar(i),~,~,Dens] = Demography(New,ParTree, ParLiana,'direct', x,mod);
            B(i)=Dens.B;        % extract the biomass for plotting
        end

    % make subplots
    subplot(422) 
        % conversion from x to D dimension for plotting
        gD = gL/((Dem.theta-1).*1.6^(ParLiana.c*(Dem.theta-1)-1)...
                *ParLiana.beta^(Dem.theta-1)*ParLiana.c);
        plot(Qrange*365, gD,LTY(p),'color', col.L, 'LineWidth',2); hold on
        set(gca,'xticklabel',[]); 
        ylabel({'Diam. growth rate';'(cm yr^-^1)'}); 
        title('Drought Shift')
        set(gca,'YAxisLocation','right')
        ylim([0.047 0.055]); 

subplot(424)
        % conversion from x to D dimension for plotting
        gD = gT/((Dem.theta-1).*33^(ParTree.c*(Dem.theta-1)-1)...
                *ParTree.beta^(Dem.theta-1)*ParTree.c);
        
        plot(Qrange*365, gD,LTY(p),'color', col.T, 'LineWidth',2); hold on
        set(gca,'xticklabel',[]); 
        ylabel({'Diam. growth rate';'(cm yr^-^1)'});  
        set(gca,'YAxisLocation','right');
        ylim([0.47 0.54]);

 subplot(426)
        plot(Qrange*365, PL*100, LTY(p), 'LineWidth',2, 'color','k'); hold on;              
        set(gca,'xticklabel',[]); 
        ylabel({'Liana prevalence';' (P_L, %)'});
        set(gca,'YAxisLocation','right')
        ylim([17 27]); 

subplot(428)
        plot(Qrange*365, B,LTY(p), 'LineWidth',2, 'color','k'); hold on
        xlabel('Dry season length');         
        ylabel({'Biomass';'(Mg ha^{-1})'})
        set(gca,'YAxisLocation','right')
        ylim([240 360]); 
        pause(.1)
end        

% add sublabels
sublabel('FontName','Helvetica','FontWeight','bold', 'BackgroundColor','none')

% add a tabular legend
    subplot(428)
    labs=[sprintf('%1.1f',beta(1));sprintf('%1.1f',beta(2)); sprintf('%1.1f',beta(3))];
    labels1 = {labs(1,:)' labs(2,:)' labs(3,:)};
    labels2 = {labs(1,:)' labs(3,:)' labs(2,:)};
    labels = cellfun(@(x,y) {sprintf('%5s%5s', x, y)}, labels1, labels2);
    l = legend(labels, 'FontName', 'FixedWidth');%, 
    titles = {''; '\beta_{T}'; '\beta_{L}'};
    titles = sprintf('%8s%8s%12s', titles{:});
    l.Title.String = titles;
    legend('boxoff')  



%% Figure 5
% Transient dynamic analyses showing how a gradual 50-year step-increase 
% in CO2 affects liana prevalence and forest structure.

% NOTE: this run will take a bit of time to find stability in the spinup!
clear;   clc; 
run 'TransientCO2.m'

% in this function we follow the following steps
% 1.    Define the CO2 scenario we want to study
% 2.    Define the distinct carbon aquisition rates under these distinct CO2
%       conditions
% 3.    Translate carbon aquisition in corresponding growth rates
% 4.    Benchmark for the starting conditions (i.e., in 1975) 
% 5.    Run the model 
%           run spinup followed by increase and then steady state for
%           recovery
% 6.    Plotting of the figure
% 7.    Defining the functions of the transient model

%% Figure 6
% Model benchmarking to physiological trait data for lianas and trees 
% measured during wet and dry season by Smith-Martin et al., 
% (2019, Parque Municipal Summit, Panama).
clear;   clc; 
run 'Analysis_Benchmarking.m'

% this script will do the following
% 1. Load the data of Smith-Martin et al, 2019 (access via Dryad)
% 2. Optimize the Physiological model to provide the best fit to KMAX, 
%       VCmax and PsiLeaf of Smith Martin et al.
% 3. Simulating the Physiological model with these optimal values
% 4. Plot the results for visual validation



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SUPPLEMENTARY MATERIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we can provide all the processing code for the analysis discussed in the
% supplementary, but due to author restrictions, we can only provide a
% functional code upon written proof of access to the data.


