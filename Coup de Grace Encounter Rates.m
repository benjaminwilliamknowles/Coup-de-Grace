clear; close all

% plot specifications
fsize=36;
labs={'calcified cells';'naked cells'};

% size and density estimates for calcified cells, naked cells, and liths
dh=[6 5 ]./10^4; % diameter of host cells (cm)
dv=.2/10^4;         % diameter of virus (cm)
rh=dh./2;           % host cell radius (cm)
rv=dv/2;            % virus radius  (cm)
re=rh+rv;           % encounter radius (cm)
rhop=[1.19 1.05 ];   % particle density (g/cm^3) from Paasche et al. 2002

% physical properties 
rhof=1.027;                 % fluid (seawater) density (g/cm^3)
nu=0.01;                    % kinematic viscosity (cm^2/s)
g=980;                      % gravitational acceleration (cm/s^2)
eta=rhof*nu;                % dynamic viscosity [g/cm/s] = rhof*kinematic viscosity [cm^2/s]

% TKE dissipation rate (cm^2/s^3)
Es=10^-4; % "still" or "calm"
Et=10^-0; % "turbulent" or "stormy"

% sinking velocities 
wh = 2/9*g*rh.^2.*(rhop-rhof)./eta; %cell velocity by Stokes' law (cm/s)
whperd=wh/100*24*3600;              %cell velocity by Stokes' law (m/d)

% concentrations
vconc=10.^[3:.1:6]'; % virus concentration (#/cm^3)
hconc=10.^[1:.1:6]'; % host  concentration (#/cm^3)
% Vconc=repmat(vconc,1,length(hconc));
% hconc=repmat(hconc,length(vconc),1);

% estimate Brownian motion kernel beta_br (cm^3/s)
k=1.38064852*10^-23;    % Boltzmann's constant (m^2 kg s^-2 K)
kcm=k*10^4*10^3;        % Boltzmann's constant (cm^2 g s^-2 K)
T=293.15;               % absolute temperature (Kelvin) at 20 oC
betabr=2.*kcm.*T.*re.^2./3./eta./rh./rv; % Brownian motion kernel (cm^3/s)

% estimate encounter kernels due to sinking and turbulence (cm^3/s)
for i=1:length(dh);
    betab(i,1)=pi.*re(i).^2.*wh(i); %encounter kernel due to sinking

    betat(i,:)=.42*pi.*re(i).^3.*sqrt(Et./nu); %encounter kernel due to strong turbulence
    betas(i,:)=.42*pi.*re(i).^3.*sqrt(Es./nu); %encounter kernel due to weak turbulence
end

%convert to volume encountered per day
betab=betab*24*3600;
betat=betat*24*3600;
betas=betas*24*3600;
betabr=betabr*24*3600;

%%%%%%%%%%%%%%%%%%%%%%%%%
bstill=betabr+betab+betas;
bturb =betabr+betab+betat;

for i=1:length(dh);
    encS(:,i)=bstill(i).*hconc;
    encT(:,i)=bturb(i).*hconc;
end

% 
% figure; hold on
% subplot(1,2,2);
% loglog(hconc,encS(:,1),hconc,encT(:,1),'linewidth',6);hold on; 
% title('calcified','fontsize',fsize)
% 
% h=legend('still','turbulence','Location','NorthWest','AutoUpdate','off')
% legend('boxoff')
% 
% set(gca,'fontsize',fsize,'xtick',10.^[1:6],'ytick',10.^[-5:2:1],'plotboxaspectratio',[1.2 1 1])
% 
% xlabel('[{\it{Ehux}}] (# ml^{-1})');
% 
% axis(10.^[ 1 6 -5 2])
% 
% vline(10^5,'k--')
% 
% subplot(1,2,1);
% loglog(hconc,encS(:,2),hconc,encT(:,2),'linewidth',6);hold on; 
% title('naked','fontsize',fsize)
% 
% set(gca,'fontsize',fsize,'xtick',10.^[1:6],'ytick',10.^[-5:2:1],'plotboxaspectratio',[1.2 1 1])
% 
% axis(10.^[1 6 -5 2 ])
% xlabel('[{\it{Ehux}}] (# ml^{-1})');
% ylabel('encounters virus^{-1} d^{-1}')
% 
% 
% vline(10^5,'k--')

% export_fig encounterspervirusperday -png -m5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on
subplot(1,2,2);
loglog(hconc,1./encS(:,1),hconc,1./encT(:,1),'linewidth',6);hold on; 
title('calcified','fontsize',fsize)

h=legend('still','turbulence','Location','SouthWest','AutoUpdate','off')
legend('boxoff')

set(gca,'fontsize',fsize,'xtick',10.^[1:6],'ytick',10.^[-2:2:5],'plotboxaspectratio',[1.2 1 1])

xlabel('[{\it{Ehux}}] (# ml^{-1})');

axis(10.^[ 1 6 -2 5])

plot(10.^[1 6],[3 3],'k--','linewidth',1.5)

% vline(10^5,'k--')

subplot(1,2,1);
loglog(hconc,1./encS(:,2),hconc,1./encT(:,2),'linewidth',6);hold on; 
title('naked','fontsize',fsize)

set(gca,'fontsize',fsize,'xtick',10.^[1:6],'ytick',10.^[-2:2:5],'plotboxaspectratio',[1.2 1 1])

axis(10.^[1 6 -2 5])
xlabel('[{\it{Ehux}}] (# ml^{-1})');
ylabel('time to encounter a host (d)')

plot(10.^[1 6],[3 3],'k--','linewidth',1.5)

% vline(10^5,'k--')

export_fig timetofindahost -png -m5


