%% Obtaining and processing of environmental data of BCI
dat = load('C:\Users\hanne\Dropbox\Colab_Matteo\Hydromodel\data\BCI_EC_v50.mat', 'jday', 'vpd' ,  'tair','swc', 'Par_tot', 'dry','tdr' );


% calculate soil water potential
            % estimated from Kupers (2019)
            a = -log(20)/0.27273;
            b = log(0.07)-a;
            n = length(dat.jday);
            psi=nan(n,4);

            for i=1:3
            psi(:,i) = -exp(a*dat.tdr(:,i)/max(dat.tdr(:,i))+b)+0.05;
            end
            psi(:,4) =nanmean(psi(:,3));

            beta = 0.961;f=zeros(4,1); %Jackson 1996
            f(1) = 1-exp(log(beta)*10);
            f(2) = exp(log(beta)*10)-exp(log(beta)*40);
            f(3) = exp(log(beta)*40)-exp(log(beta)*100);
            f(4) = exp(log(beta)*100);

            %plot(psi*f)
            dat.psiS = psi*f;

%% plot figure
figure(1); clf;
pt=2;
tod=linspace(0,24-5/60, 288);

subplot(221)            % plot PsiS
        dry=dat.psiS(dat.dry);
        wet=dat.psiS(~dat.dry);
                PsiS(1)=nanmedian(dry);
                bar(1, PsiS(1), 'red');        hold on;
                 iqt = quantile(dry,[0.25 0.75]) - PsiS(1);
                     er = errorbar(1,PsiS(1), iqt(1),iqt(2));    
                     er.Color = [0 0 0];                            
                     er.LineStyle = 'none'; 
              
                PsiS(2)=nanmedian(wet);
                bar(2,  PsiS(2), 'blue');        hold on;
                 iqt = quantile(wet,[0.25 0.75]) - PsiS(2);
                     er = errorbar(2,PsiS(2), iqt(1),iqt(2));    
                     er.Color = [0 0 0];                            
                     er.LineStyle = 'none'; 
                     
         ylabel(['\Psi_S']);
         set(gca,'XTick',[]);
         set(gca,'XTick',[]);
         text(1,.1,'DRY','HorizontalAlignment','center','FontWeight', 'Bold');
         text(2,.1,'WET','HorizontalAlignment','center','FontWeight', 'Bold');
         
        subplot(222)
        % reshape the data
        D = zeros(288,2);               % VPD
        [~,D(:,1)] =mdv(dat.vpd,288,dat.dry);
        [~,D(:,2)]=mdv(dat.vpd,288, ~dat.dry);
        plot(tod,D(:,1),'color','red','LineWidth',2); hold on;
        plot(tod,D(:,2),'color','blue','LineWidth',2); 
        ylabel('D (kPa)')
       
        
        subplot(223)
        I0 = zeros(288,2);              % Incoming radiation
        [~,I0(:,1)] =mdv(dat.Par_tot,288,dat.dry);
        [~,I0(:,2)]=mdv(dat.Par_tot,288, ~dat.dry);
        plot(tod,I0(:,1),'color','red','LineWidth',2); hold on;
        plot(tod,I0(:,2),'color','blue','LineWidth',2); 
        ylabel('I_0 (\mumol m^-^2)');       xlabel('TOD (h)');

        
        subplot(224)
        Tc = zeros(288,2);          % temperature at leaf level
        [~,Tc(:,1)] =mdv(dat.tair, 288, dat.dry);
        [~,Tc(:,2)]=mdv(dat.tair, 288, ~dat.dry);
        plot(tod,Tc(:,1),'color','red','LineWidth',2); hold on;
        plot(tod,Tc(:,2),'color','blue','LineWidth',2); 
        ylabel('T_c (C)');       xlabel('TOD (h)');

savefig('AvgDay.fig')


%% write it all away
% the second column is always the wet conditions
DAY.Time = tod;
DAY.W.Tc=Tc(:,2);
DAY.D.Tc=Tc(:,1);
DAY.W.I0=I0(:,2);
DAY.D.I0=I0(:,1);
DAY.W.D=D(:,2);
DAY.D.D=D(:,1);
DAY.W.PsiS=repelem(PsiS(:,2), length(tod));
DAY.D.PsiS=repelem(PsiS(:,1), length(tod));
DAY.DSL=sum(dat.dry)/length(dat.dry);

save('AVGday.mat','DAY')