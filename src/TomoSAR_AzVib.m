% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
clear;clc;%close all
addpath('./SES_functions')
c = physconst('LightSpeed');
export_flag = 0;
export_directory = '';
if export_flag
    mkdir(export_directory)
    export_directory=[export_directory,'Sim-AzVib\',datestr(now,'yyyymmdd_HHMM'),'\'];
    mkdir(export_directory);
end
fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);
%% Parameters:
fc = 77e9; lambda = c/fc;
bw = 1e9;
T = 60e-6; cr = bw/T;
prt = 0.08; % (s)
snr = -10;25;
Nr = 256;512;
Na = 8;
Ne = 1;
epochs = 50;200;500;
d_az = lambda/2;
d_el = lambda/2;
time_tot = prt*(epochs-1); % Total monitoring time

signal_model = 'polar'; % 'cartesian' or 'polar'
mode_az='mimo';
mode_el='mimo';
if strcmp('CARTESIAN',upper(signal_model))
    mode_el=mode_az; % Use similar imaging modes when using 'cartesian' signal model!!
end
if strcmp('MIMO',upper(mode_az))
    m_az=1;
elseif strcmp('MONO',upper(mode_az))
    m_az=2;
end
if strcmp('MIMO',upper(mode_el))
    m_el=1;
elseif strcmp('MONO',upper(mode_el))
    m_el=2;
end

%% Radar position:
X_rad=d_az*(0:Na-1);
Y_rad=0;
Z_rad=0*(1:epochs);

%% Target position:
% Polar
R_tar = 1*[20,20,20]; theta_tar = [-10,0,20]; ph_tar=zeros(size(R_tar));
% R_tar = .5*[20]; theta_tar = 0;[5.2]; ph_tar=zeros(size(R_tar));
% Cartesian:
X_tar=R_tar.*sind(theta_tar);Y_tar=R_tar.*sqrt(1-sind(theta_tar).^2-sind(ph_tar).^2);Z_tar=R_tar.*sind(ph_tar);
num_tar=length(R_tar);
 
%% Displacement functions:
% Fluctuations:
dR_amp_tar = [0.005,0.001,0.002];%,.00015];
% dR_amp_tar = 10e-4;%2.324e-4;%,.00015];
dR_f_tar = [7,2,.4];%./prt;%,0.01]; % displacement frequency Hz(1/epochs)
% dR_f_tar = 5;%24.5;%24.76855;%./prt;%,0.01]; % displacement frequency Hz(1/epochs)
vib_mode = 'SINE'; % Vibration signal model: 'sine' or 'triangle'
if strcmp('SINE',upper(vib_mode))
    dR_tar_fluc = displacement_model_sin_TS( dR_amp_tar,dR_f_tar,prt,time_tot );
elseif strcmp('TRIANGLE',upper(vib_mode))
    dR_tar_fluc = displacement_model_triangle_TS( dR_amp_tar,dR_f_tar,prt,time_tot );
end
% Total displacement:
dR_tar_total_disp = dR_tar_fluc+(5e-4)*(rand(epochs,1)-0.5);
figure;plot(dR_tar_total_disp,'LineWidth',1.5);xlabel('Time samples');ylabel('Displacement (m)')

%% SIGNAL
TomoSAR_SIGNAL

%% Show parameters in a Table
table_row = {'fc (GHz)';'Wavelength (m)';'Bandwidth (GHz)';'Chirp duration (s)';'snr (dB)';'Range samples (Nr)';'Azimuth samples (Na)';'Azimuth antenna spacing (m)';'Azimuth mode';...
             'Azimuth beam (deg)';'epochs';'Targets number';'Targets Range (m)';'Targets Phi (deg)';'Targets Theta (deg)';'Targets Vib. Amp (m)';'Targets Vib. f (1/ep)'};
table_val = {fc/1e9;lambda;bw/1e9;T;snr;Nr;Na;d_az;upper(mode_az);beam_az;epochs;num_tar;num2str(R_tar);num2str(ph_tar);num2str(theta_tar);num2str(dR_amp_tar);num2str(dR_f_tar)};
Table = table(table_val,'RowNames',table_row);
figure;tbl=uitable('Data', table_val, 'RowName', table_row,'ColumnName',{'Value'},'Position',[10,10,500,400]);tbl.ColumnWidth={num_tar*50};
if export_flag
    print(gcf,[export_directory , 'Table of Parameters.jpg'],'-djpeg','-r400');
end
%% ========================================================================
%% +-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*
%% +-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*
%% ========================================================================

%% Range compression:
rc=fft(signal_MIMO_TS,[],1);
slc =fftshift(fft(rc,[181],2),2);
[r_tar,c_tar] = find( abs(slc(:,:,1))==max(max(abs(slc(:,:,1)))) );
[~,~,TS_cum_phase] = TSInSAR(slc);
defo_map_TS = phase2defo(TS_cum_phase,lambda);
dif1=squeeze(defo_map_TS(r_tar,c_tar,:));
time_ax=prt*(1:1:length(dif1));
figure;
subplot(2,1,1);plot(time_ax(:),dif1(:));xlabel('time')
f_dif1 = fft(dif1-mean(dif1));
r_freq = 1/(prt*length(f_dif1)); % frequency resolution!
f_ax = r_freq*(0:1:length(f_dif1)-1);
subplot(2,1,2);plot(f_ax,abs(f_dif1));xlabel('f')

% slc_vel =fftshift(fft(slc,[146],3),3);
% figure;subplot(1,2,1);imagesc(abs(slc(:,:,1)));xlabel('Az');ylabel('R');
% subplot(1,2,2);imagesc(abs(squeeze(slc_vel(r_tar,:,:))));xlabel('Vel');ylabel('Az');


%% Spectral estimation:
w_r = [r_tar-0,r_tar+0];
time_vec = (0:epochs-1).*prt;
antenna_loc = d_az*(0:Na-1)';
% theta = 0;-0:2:10;%-round(beam_az/2):1:round(beam_az/2); 
% fdisp = [0:.5:30];(0:.025:.25)./prt; % displacement frequency
% Adisp = (0:5e-5:20e-4);(0:0.0005:0.008); % displacement amplitude
theta = [-30:2:30]; %round(theta_tar)*[0.5:0.1:1.5]; 
fdisp = [ 0:.2:12 ]; % displacement frequency
Adisp = linspace( 0,8,17)*1e-3; % displacement amplitude
%
% [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVib(rc,w_r,lambda,time_vec,antenna_loc,theta,Adisp,fdisp,num_tar,mode_az,vib_mode);
c_iter = 4;
[Pbf,Pcapon,Pmusic,Pbf_C,Pcapon_C,Pmusic_C]=SpectralEstimation_AzVib_CLEAN(rc,w_r,lambda,time_vec,antenna_loc,theta,Adisp,fdisp,num_tar,mode_az,vib_mode,c_iter);
%

% Pbf = absdB(Pbf); Pcapon = absdB(Pcapon); Pmusic = absdB(Pmusic);
Pbf = norm01((Pbf)); Pcapon = norm01((Pcapon)); Pmusic = norm01((Pmusic));
Pbf_C = norm01((Pbf_C)); Pcapon_C = norm01((Pcapon_C)); Pmusic_C = norm01((Pmusic_C));
Pbf = absdB(Pbf); Pcapon = absdB(Pcapon); Pmusic = absdB(Pmusic);
Pbf_C = absdB(Pbf_C); Pcapon_C = absdB(Pcapon_C); Pmusic_C = absdB(Pmusic_C);
caxislim=[-10,0]; 
[rrr,cc,vv]=ind2sub(size(Pbf),find(Pbf==1));       thet_bf    = theta(rrr),Ampd_bf   =Adisp(cc), freq_bf    = fdisp(vv)
[rrr,cc,vv]=ind2sub(size(Pcapon),find(Pcapon==1)); thet_capon = theta(rrr),Ampd_capon=Adisp(cc), freq_capon = fdisp(vv)
[rrr,cc,vv]=ind2sub(size(Pmusic),find(Pmusic==1)); thet_music = theta(rrr),Ampd_music=Adisp(cc), freq_music = fdisp(vv)


snr_bf    = 10*log10( 1/std(Pbf(15:end,15:end,8:end),[],'all') );
snr_capon = 10*log10( 1/std(Pcapon(15:end,15:end,8:end),[],'all') );
snr_music = 10*log10( 1/std(Pmusic(15:end,15:end,8:end),[],'all') );

for tar_qq = 1:num_tar
    thet_idx (tar_qq) = find(round(theta) == round(theta_tar (tar_qq)));
    Adisp_idx(tar_qq) = find(round(Adisp,4) == round(dR_amp_tar(tar_qq),4));
    fdisp_idx(tar_qq) = find(round(fdisp,1) == round(dR_f_tar  (tar_qq),1));
end
P_idx = sub2ind([length(theta), length(Adisp), length(fdisp)],thet_idx,Adisp_idx,fdisp_idx);
snr_bf      = 10*log10( mean(Pbf     (P_idx))/((sum( Pbf     (:) )-sum( Pbf     ( P_idx ) ))/(length(Pbf     (:))-num_tar)) );
snr_capon   = 10*log10( mean(Pcapon  (P_idx))/((sum( Pcapon  (:) )-sum( Pcapon  ( P_idx ) ))/(length(Pcapon  (:))-num_tar)) );
snr_music   = 10*log10( mean(Pmusic  (P_idx))/((sum( Pmusic  (:) )-sum( Pmusic  ( P_idx ) ))/(length(Pmusic  (:))-num_tar)) );
snr_bf_C    = 10*log10( mean(Pbf_C   (P_idx))/((sum( Pbf_C   (:) )-sum( Pbf_C   ( P_idx ) ))/(length(Pbf_C   (:))-num_tar)) );
snr_capon_C = 10*log10( mean(Pcapon_C(P_idx))/((sum( Pcapon_C(:) )-sum( Pcapon_C( P_idx ) ))/(length(Pcapon_C(:))-num_tar)) );
snr_music_C = 10*log10( mean(Pmusic_C(P_idx))/((sum( Pmusic_C(:) )-sum( Pmusic_C( P_idx ) ))/(length(Pmusic_C(:))-num_tar)) );

set(0,'DefaultAxesFontName','times','DefaultTextFontName','times','defaultAxesFontSize',12);
figure('Position', [30 60 1500 700]); %sgtitle('Azimuth/Vibration');
subplot(3,3,1);imagesc(fdisp,1e3*Adisp,squeeze(max(Pbf,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['BF | SNR: ',num2str(snr_bf)],'','\xi (mm)'});  xlabel('\eta (Hz)');title('\eta - \xi')
subplot(3,3,2);imagesc(fdisp,theta,squeeze(max(Pbf,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);    ylabel('\theta (deg)');         xlabel('\eta (Hz)');title('\eta - \theta')
subplot(3,3,3);imagesc(1e3*Adisp,theta,max(Pbf,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');         xlabel('\xi (mm)');title('\xi - \theta')

subplot(3,3,4);imagesc(fdisp,1e3*Adisp,squeeze(max(Pcapon,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['Capon | SNR: ',num2str(snr_capon)],'','\xi (mm)'}); xlabel('\eta (Hz)');
subplot(3,3,5);imagesc(fdisp,theta,squeeze(max(Pcapon,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('\eta (Hz)');
subplot(3,3,6);imagesc(1e3*Adisp,theta,max(Pcapon,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('\xi (mm)');

subplot(3,3,7);imagesc(fdisp,1e3*Adisp,squeeze(max(Pmusic,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);xlabel('\eta (Hz)');ylabel({['MUSIC | SNR: ',num2str(snr_music)],'','\xi (mm)'})
subplot(3,3,8);imagesc(fdisp,theta,squeeze(max(Pmusic,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);xlabel('\eta (Hz)');ylabel('\theta (deg)')
subplot(3,3,9);imagesc(1e3*Adisp,theta,max(Pmusic,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);xlabel('\xi (mm)'); ylabel('\theta (deg)')
cbar = colorbar('Position',[0.93 0.11 0.02 0.81]); cbar.Title.String = "[dB]";
Lgnd = legend({'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.9; Lgnd.Position(2) = 0.01;

if export_flag
    print(gcf,[export_directory , '2D normalized PSD.jpg'],'-djpeg','-r400');
end
% CLEAN:
set(0,'DefaultAxesFontName','times','DefaultTextFontName','times','defaultAxesFontSize',12);
figure('Position', [30 60 1500 700]); %sgtitle('Azimuth/Vibration');
subplot(3,3,1);imagesc(fdisp,1e3*Adisp,squeeze(max(Pbf_C,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['BF | SNR: ',num2str(snr_bf_C)],'','\xi (mm)'});  xlabel('\eta (Hz)');title('\eta - \xi')
subplot(3,3,2);imagesc(fdisp,theta,squeeze(max(Pbf_C,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);    ylabel('\theta (deg)');         xlabel('\eta (Hz)');title('\eta - \theta')
subplot(3,3,3);imagesc(1e3*Adisp,theta,max(Pbf_C,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');         xlabel('\xi (mm)');title('\xi - \theta')

subplot(3,3,4);imagesc(fdisp,1e3*Adisp,squeeze(max(Pcapon_C,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['Capon | SNR: ',num2str(snr_capon_C)],'','\xi (mm)'}); xlabel('\eta (Hz)');
subplot(3,3,5);imagesc(fdisp,theta,squeeze(max(Pcapon_C,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('\eta (Hz)');
subplot(3,3,6);imagesc(1e3*Adisp,theta,max(Pcapon_C,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('\xi (mm)');

subplot(3,3,7);imagesc(fdisp,1e3*Adisp,squeeze(max(Pmusic_C,[],1)));caxis(caxislim);hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);xlabel('\eta (Hz)');ylabel({['MUSIC | SNR: ',num2str(snr_music_C)],'','\xi (mm)'})
subplot(3,3,8);imagesc(fdisp,theta,squeeze(max(Pmusic_C,[],2)));caxis(caxislim);    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);xlabel('\eta (Hz)');ylabel('\theta (deg)')
subplot(3,3,9);imagesc(1e3*Adisp,theta,max(Pmusic_C,[],3));caxis(caxislim);         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);xlabel('\xi (mm)'); ylabel('\theta (deg)')
cbar = colorbar('Position',[0.93 0.11 0.02 0.81]); cbar.Title.String = "[dB]";
Lgnd = legend({'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.9; Lgnd.Position(2) = 0.01;

if export_flag
    print(gcf,[export_directory , '2D normalized PSD_CLEAN.jpg'],'-djpeg','-r400');
end





% set(0,'DefaultAxesFontName','times','DefaultTextFontName','times','defaultAxesFontSize',12);
% figure('Position', [30 60 1500 700]); sgtitle('Azimuth/Vibration');
% subplot(4,3,1);imagesc(fdisp,1e3*Adisp,squeeze(max(Pmf,[],1)));hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['MF | SNR: ',num2str(snr_mf)],'','Vib. A (m)'}); xlabel('Vib. f (Hz)');title('Vib. f - Vib. A')
% subplot(4,3,2);imagesc(fdisp,theta,squeeze(max(Pmf,[],2)));    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);    ylabel('\theta (deg)');         xlabel('Vib. f (Hz)');title('Vib. f - \theta')
% subplot(4,3,3);imagesc(1e3*Adisp,theta,max(Pmf,[],3));         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');         xlabel('Vib. A (mm)');title('Vib. A - \theta')
% cbar = colorbar('Position',[0.93 0.69 0.02 0.2]);
% 
% subplot(4,3,4);imagesc(fdisp,1e3*Adisp,squeeze(max(Pbf,[],1)));hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['BF | SNR: ',num2str(snr_bf)],'','Vib. A (m)'}); xlabel('Vib. f (Hz)');title('Vib. f - Vib. A')
% subplot(4,3,5);imagesc(fdisp,theta,squeeze(max(Pbf,[],2)));    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);    ylabel('\theta (deg)');         xlabel('Vib. f (Hz)');title('Vib. f - \theta')
% subplot(4,3,6);imagesc(1e3*Adisp,theta,max(Pbf,[],3));         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');         xlabel('Vib. A (mm)');title('Vib. A - \theta')
% cbar = colorbar('Position',[0.93 0.69 0.02 0.2]);
% 
% subplot(4,3,7);imagesc(fdisp,1e3*Adisp,squeeze(max(Pcapon,[],1)));hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);ylabel({['Capon | SNR: ',num2str(snr_capon)],'','Vib. A (mm)'});xlabel('Vib. f (Hz)');
% subplot(4,3,8);imagesc(fdisp,theta,squeeze(max(Pcapon,[],2)));    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('Vib. f (Hz)');
% subplot(4,3,9);imagesc(1e3*Adisp,theta,max(Pcapon,[],3));         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel('Vib. A (mm)');
% cbar = colorbar('Position',[0.93 0.4 0.02 0.2]);
% 
% subplot(4,3,10);imagesc(fdisp,1e3*Adisp,squeeze(max(Pmusic,[],1)));hold on;scatter(dR_f_tar',1e3*dR_amp_tar',150, 'ro','LineWidth',1.5);xlabel('Vib. f (Hz)');ylabel({['MUSIC | SNR: ',num2str(snr_music)],'','Vib. A (mm)'})
% subplot(4,3,11);imagesc(fdisp,theta,squeeze(max(Pmusic,[],2)));    hold on;scatter(dR_f_tar',theta_tar',150,  'ro','LineWidth',1.5);xlabel('Vib. f (Hz)');ylabel('\theta (deg)')
% subplot(4,3,12);imagesc(1e3*Adisp,theta,max(Pmusic,[],3));         hold on;scatter(1e3*dR_amp_tar',theta_tar',150,'ro','LineWidth',1.5);xlabel('Vib. A (mm)'); ylabel('\theta (deg)')
% cbar = colorbar('Position',[0.93 0.11 0.02 0.2]);%cbar.Title.String = "[dB]";
% Lgnd = legend({'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.9; Lgnd.Position(2) = 0.01;

if export_flag
    print(gcf,[export_directory , '2D normalized power spectrum density.jpg'],'-djpeg','-r400');
end
% % Point cloud:
figure('Position', [30 60 1500 500]); sgtitle('Point Cloud')
th_mf = .8; % Threshold
pc_mf = gen_pointcloud(Pmf,th_mf);thet = theta(pc_mf(:,1));Ampd=1e3*Adisp(pc_mf(:,2));freq = fdisp(pc_mf(:,3)); power = pc_mf(:,4);
subplot(2,2,1);scatter3(thet,Ampd,freq,[],power,'filled'); hold on; scatter3(theta_tar',1e3*dR_amp_tar',dR_f_tar',200,'r');
colormap(jet);caxis([th_mf,1]);title('BF'); zlabel('Vib. f (Hz)');xlabel('\theta (deg)');ylabel('Vib. A (mm)');zlabel('Vib. f (Hz)');
xlim([min(theta),max(theta)]);ylim(1e3*[min(Adisp),max(Adisp)]);zlim([min(fdisp),max(fdisp)])

th_bf = .8; % Threshold
pc_bf = gen_pointcloud(Pbf,th_bf);thet = theta(pc_bf(:,1));Ampd=1e3*Adisp(pc_bf(:,2));freq = fdisp(pc_bf(:,3)); power = pc_bf(:,4);
subplot(2,2,2);scatter3(thet,Ampd,freq,[],power,'filled'); hold on; scatter3(theta_tar',1e3*dR_amp_tar',dR_f_tar',200,'r');
colormap(jet);caxis([th_bf,1]);title('BF'); zlabel('Vib. f (Hz)');xlabel('\theta (deg)');ylabel('Vib. A (mm)');zlabel('Vib. f (Hz)');
xlim([min(theta),max(theta)]);ylim(1e3*[min(Adisp),max(Adisp)]);zlim([min(fdisp),max(fdisp)])

th_capon = .80; % Threshold
pc_capon = gen_pointcloud(Pcapon,th_capon);thet = theta(pc_capon(:,1));Ampd=1e3*Adisp(pc_capon(:,2));freq = fdisp(pc_capon(:,3)); power = pc_capon(:,4);
subplot(2,2,3);scatter3(thet,Ampd,freq,[],power,'filled'); hold on; scatter3(theta_tar',1e3*dR_amp_tar',dR_f_tar',200,'r');
colormap(jet);caxis([th_capon,1]);title('Capon');zlabel('Vib. f (Hz)');xlabel('\theta (deg)');ylabel('Vib. A (mm)');zlabel('Vib. f (Hz)');
xlim([min(theta),max(theta)]);ylim(1e3*[min(Adisp),max(Adisp)]);zlim([min(fdisp),max(fdisp)])

th_music = .80; % Threshold
pc_music = gen_pointcloud(Pmusic,th_music);thet = theta(pc_music(:,1));Ampd=1e3*Adisp(pc_music(:,2));freq = fdisp(pc_music(:,3)); power = pc_music(:,4);
subplot(2,2,4);scatter3(thet,Ampd,freq,[],power,'filled'); hold on; gt=scatter3(theta_tar',1e3*dR_amp_tar',dR_f_tar',200,'r');
colormap(jet);cbar = colorbar();caxis([th_music,1]);title('MUSIC');xlabel('\theta (deg)');ylabel('Vib. A (mm)');zlabel('Vib. f (Hz)');
xlim([min(theta),max(theta)]);ylim(1e3*[min(Adisp),max(Adisp)]);zlim([min(fdisp),max(fdisp)])
cbar = colorbar('Position',[0.93 0.168 0.022 0.7]);%cbar.Title.String = "dB";
colormap(cbar,'jet')
Lgnd = legend(gt,{'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.9; Lgnd.Position(2) = 0.05;
if export_flag
    print(gcf,[export_directory , '3D Point Cloud.jpg'],'-djpeg','-r400');
end
