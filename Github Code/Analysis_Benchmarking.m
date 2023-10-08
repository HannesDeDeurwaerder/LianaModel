%% VALIDATION of model performance
load('ModelParameterV1.0.mat') 
load('AVGday.mat')  

% declare constants
rho =     1e3;               % density of water [kg m-3]                            
g =       9.81;              % gravity [m s-2]
cf=2*12e-9*3600*24*365;      % conversion factor from umol CO2 m-2 s-1 to Kg dry mass yr-1
labels={'LW','TW','LD','TD'};% labels: liana wet, tree wet, liana dry, & tree dry

% declare matrices
An=nan(1,4);               Psi_md=nan(1,4);
Psi_pd=nan(1,4);           Psi_diff=nan(1,4);
Gs=nan(1,4);

%% LOAD AND PREPARE DATASET OF SMITH-MARTIN ET AL. 2019
% data is accessible via Dryad database

% load validation dataset
DS=readtable('.\SmithMartin_2019.csv');

    DS(isnan(DS.gs),:) = [];            % remove values that did not measure gs and An on same leaves
    DS(isnan(DS.Psi_pd),:) = [];        % remove lines that have no Predawn
    
    DS.Psi_diff=DS.Psi_md-DS.Psi_pd;    % add column with Psi Diff
    DS.dry=DS.month<5;                  % indicate the dry season
    habit=char(DS.Habit);    DS.liana=habit(:,1)=='l';      % indicate the lianas
    
    % groups needed for the boxplot analysis
    DS.group(habit(:,1)=='l' & ~DS.dry)=1;    % Liana wet
    DS.group(habit(:,1)=='l' & DS.dry)=2;     % Liana dry
    DS.group(habit(:,1)=='t' & ~DS.dry)=3;    % Tree wet
    DS.group(habit(:,1)=='t' & DS.dry)=4;     % Tree dry

    
    % calculate the predawn water potentials at Smith martin et al location
    Psi_pd(1)=mean(DS.Psi_pd(DS.liana & ~DS.dry )) ; % liana wet
    Psi_pd(2)=mean(DS.Psi_pd(DS.liana & DS.dry ));   % liana dry
    Psi_pd(3)=mean(DS.Psi_pd(~DS.liana & ~DS.dry )); % tree wet
    Psi_pd(4)=mean(DS.Psi_pd(~DS.liana & DS.dry ));  % tree dry

    % calculate the soil water potentials at smith martin et al location
    % NOTE that a gravimetric term should be accounted
    % we estimate forest height of sampling assuring a minimal 
    % negative psi soil during wet conditions
    height = 1;                            % estimated height
    PsiSoil = Psi_pd;

    % used temperature and radiance in the gass exchange
    temp=32;    I0=1200;
    % assumed VPD conditions  in dry and wet season, based on observations in BCI
    vpdd=1.7;   vpdw=0.6;

 %% OPTIMIZER TO KMAX, VCMAX AND PSILEAF
        % ranges of the values being optimized on
        ymin=[1 10 0.5];    ymax=[15 100 10];  
                % 1st is Kmax
                % 2nd is VCmax
                % 3rd is psileaf

        % LIANA
        %------------------------------------------------------------------
        % initial estimates
        y0=[6.5   50    5]; 
        params.PAR=ParLiana;
        params.datw= struct('I',I0,'Ca',400,'TL',temp,'D',vpdw,'psiS',PsiSoil(1),'H',height,'LightLimited','n');
        params.datd= struct('I',I0,'Ca',400,'TL',temp,'D',vpdd,'psiS',PsiSoil(2),'H',height,'LightLimited','n');
        params.Psi_pd=[Psi_pd(1) Psi_pd(2)];

        % benchmark data
        use=(DS.group==1);
        params.BMwet=[median(DS.An(use)) median(DS.gs(use)) median(DS.Psi_diff(use)) ];
        use=(DS.group==2);
        params.BMdry=[median(DS.An(use)) median(DS.gs(use)) median(DS.Psi_diff(use)) ];

        % run optimizer function
        yliana = lsqnonlin(@(y) OptimLiana(y,params),y0,ymin,ymax);
        clear params
        
        % TREE
        %------------------------------------------------------------------
        % initial estimates
        y0=[6.5   40    3];     
        params.PAR=ParTree;
        params.datw= struct('I',I0,'Ca',400,'TL',temp,'D',vpdw,'psiS',PsiSoil(3),'H',height,'LightLimited','n');
        params.datd= struct('I',I0,'Ca',400,'TL',temp,'D',vpdd,'psiS',PsiSoil(4),'H',height,'LightLimited','n');
        params.Psi_pd=[Psi_pd(3) Psi_pd(4)];

        %benchmark data
        use=(DS.group==3);
        params.BMwet=[median(DS.An(use)) median(DS.gs(use)) median(DS.Psi_diff(use)) ];
        use=(DS.group==4);
        params.BMdry=[median(DS.An(use)) median(DS.gs(use)) median(DS.Psi_diff(use)) ];

        % run optimizer function
        ytree = lsqnonlin(@(y) OptimLiana(y,params),y0,ymin,ymax);
        

%% RUN THE PHYSIOLOGICAL MODEL WITH OBTAINED VALUES
%------------------------------------------------------------------------------------    

    % Assign the optimized values for the run
    ParLiana.Kmax(1:3,:)=yliana(1);         % assign the new, optimized value of Kmax for liana
    ParLiana.Vcmax25=yliana(2);             % assign the new, optimized value of VCmax for liana
    ParLiana.cost= @(x)yliana(3)*(x/ParLiana.TLP).^ParLiana.t1; % assign the new, optimized value of Psileaf for liana
    
    ParTree.Kmax(1:3,:)=ytree(1);           % assign the new, optimized value of Kmax for liana
    ParTree.Vcmax25=ytree(2);               % assign the new, optimized value of VCmax for liana
    ParTree.cost= @(x) ytree(3)*(x/ParTree.TLP).^ParTree.t1;    % assign the new, optimized value of Psileaf for liana
    
    % run the physiological model

        % wet lianas
        dat= struct('I',I0,'Ca',400,'TL',temp,'D',vpdw,'psiS',PsiSoil(1),'H',height,'LightLimited','n');
        out=An_StomatOpt(ParLiana.Vcmax25, ParLiana.Kmax, ParLiana.p50, ParLiana.gamma, dat, ParLiana.cost);
            Psi_md(1)=out.CL.psiL; 
            An(1)=out.CL.An;                
            Gs(1)=out.CL.gs;
            Psi_diff(1) = Psi_md(1)-Psi_pd(1);
        % Liana Dry
        dat= struct('I',I0,'Ca',400,'TL',temp,'D',vpdd,'psiS',PsiSoil(2),'H',height,'LightLimited','n');
        out=An_StomatOpt(ParLiana.Vcmax25, ParLiana.Kmax, ParLiana.p50, ParLiana.gamma, dat, ParLiana.cost);
            Psi_md(2)=out.CL.psiL; 
            An(2)=out.CL.An;                
            Gs(2)=out.CL.gs;
            Psi_diff(2) = Psi_md(2)-Psi_pd(2);
        % Tree Wet
        dat= struct('I',I0,'Ca',395,'TL',temp,'D',vpdw,'psiS',PsiSoil(3),'H',height,'LightLimited','n');
        out=An_StomatOpt(ParTree.Vcmax25, ParTree.Kmax, ParTree.p50, ParTree.gamma, dat, ParTree.cost);
            Psi_md(3)=out.CL.psiL; 
            An(3)=out.CL.An;                
            Gs(3)=out.CL.gs;
            Psi_diff(3) = Psi_md(3)-Psi_pd(3);           
        % Tree Dry
        dat= struct('I',I0,'Ca',395,'TL',temp,'D',vpdd,'psiS',PsiSoil(4),'H',height,'LightLimited','n');
        out=An_StomatOpt(ParTree.Vcmax25, ParTree.Kmax, ParTree.p50, ParTree.gamma, dat, ParTree.cost);
            Psi_md(4)=out.CL.psiL; 
            An(4)=out.CL.An;                
            Gs(4)=out.CL.gs;
            Psi_diff(4) = Psi_md(4)-Psi_pd(4);
            

%% PLOT FIGURE 
figure(5); clf;

    % assign to Dataset, declare graphical parameters.
    data=[DS.Psi_diff DS.An DS.gs log10(DS.An./DS.gs)];  
    ylims=[-3.2 0; 0 20; -0.2 0.7; 1 3];
    cols=['r';'b'];
    benchmark=[Psi_diff; An; Gs; log10(An./Gs)];
 
    % plot the distinct boxplots
for i =1:8
 subplot(2,4,i)
 
 if ismember(i,1:4)
     use=ismember(DS.group,[1 2]); % for uneven=only lianas
     group=DS.group(use);
     corr=0;
 else 
     use=ismember(DS.group,[3 4]); % for even = only trees
     group=DS.group(use)-2;
     corr=4;
 end 
  
 boxplot(data(use, i-corr),group,'Widths',0.6);  hold on;
 ylim(ylims(i-corr,:));
 box off;
 
% y axis visualization + provide labels
if i==1|| i==5
    set(gca,'ytick',[-3:1:0])
    
elseif i==2||i==6
    set(gca,'ytick',[0:5:20])
    
elseif i==3 || i==7
    set(gca,'ytick',[0.1:0.2:0.8])
    
elseif i==4||i==8
    set(gca,'ytick',[1:3])
    
end

set(gca, 'XColor','none')

% color the boxes
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),cols(j,:),'FaceAlpha',.1,'linestyle','none');
end

% change median line to black and thicker and change box outline
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k', 'LineWidth',2);
lines = findobj(gcf, 'type', 'line', 'Tag', 'Box');
set(lines, 'Color', 'k');

% plot model output
if ismember(i, 1:4); use=1:2; else; use=3:4; end
plot(1:2, benchmark(i-corr,use),    'MarkerFaceColor',[0.850 0.325 0.098],...
    'MarkerEdgeColor',[1 1 1], 'MarkerSize',12,...
    'Marker','pentagram',   'LineStyle','none')

% add legend
    if i==4
    h2 = get(gca,'children');
    legend(h2([2 3 1]),'dry','wet','model','fontsize',10);legend('boxoff')
    end
    
%
if i==1; ylabel('Liana','color',col.L,'fontsize',12,'fontweight','bold'); end
if i ==5; ylabel('Tree','color',col.T,'fontsize',12,'fontweight','bold'); end
end

% add the x texts
annotation(gcf,'textbox',[ 0.1953    0.8807   0.0672    0.0907], 'String',{'\Delta\Psi','(MPa)'},'linestyle','none', 'horizontalalignment','center');
annotation(gcf,'textbox',[0.3688     0.8807    0.1346    0.0907],'String',{'A_N',' (\mumol CO_2 m^-^2 s^-^1)'},'linestyle','none', 'horizontalalignment','center');
annotation(gcf,'textbox',[ 0.5632    0.8807   0.1237    0.0907], 'String',{'g_s',' (mol H_2O m^-^2 s^-^1)'},'linestyle','none', 'horizontalalignment','center');
annotation(gcf,'textbox',[  0.7492     0.8807    0.1711    0.0907], 'String',{'iWUE',' (\mumol CO_2/mol H_2O)'},'linestyle','none', 'horizontalalignment','center');

% outliers gray
lines = findobj(gcf, 'type', 'line', 'Tag', 'Outliers');
set(lines, 'MarkerEdgeColor', [1 1 1]*0.7);

% make entire figure white
set(gcf,'color',[1 1 1])

