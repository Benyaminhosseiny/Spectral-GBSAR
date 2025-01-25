% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916

% Estimating vibrational movements from single-pass SAR
%   * Latest Update: 2023-09-12 *

clear;clc;%close all
% addpath('./SES_functions')
c = physconst('LightSpeed');
export_flag = 0;
export_directory = '';
if export_flag
    mkdir(export_directory)
    export_directory=[export_directory,'Sim-AzVib-RT\',datestr(now,'yyyymmdd_HHMM'),'\'];
    mkdir(export_directory);
end
fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);

%% Parameters:
fc = 77e9; lambda = c/fc;
bw = 0.3e9;
T = 60e-6; cr = bw/T;
snr = 10;
Nr = 512;
Na = 150;100;%300;1000;
Ne = 1;
epochs = Na;
prt = 0.02;0.04;1*0.008; % (s)
v_rad = 0.025;1*0.005;
% d_az = v_rad*prt;
d_az = 1*lambda/8;
d_el = lambda/4;
time_tot = prt*(epochs-1); % Total monitoring time

signal_model = 'polar'; % 'cartesian' or 'polar'
% % USE 'cartesian' ONLY FOR EVALUATING THE IMPACTS OF MECHANICAL ERRORS! OTHERWISE USE 'polar' SIGNAL MODEL. 
mode_az='mono';
mode_el='mono';
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
X_rad=d_az*(0:Na-1); % X_rad = 1*[ 0,X_rad(randperm(length(X_rad)-1)) ];

% L_az = 80e-3; NL_az = round(L_az/d_az);
% L_stop = 1*40e-3; NL_stop = round(L_stop/d_az);
% stop_time0 = 0*2.5; % s
% stop_time = 16; % s
% 
% X_rad = zeros(1,1+round(stop_time0/prt));
% for ii = 1:floor(epochs/NL_az/2)+1
%     if stop_time*L_stop ~= 0
%         
%         for jj = 1:floor(L_az/L_stop) % increasing
%             X_rad = [X_rad,X_rad(end)+d_az*[0:NL_stop]]; % Linear
%             X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
%         end
%         
% %         X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt/2))]; %stop
%         
%         for jj = 1:floor(L_az/L_stop) % decreasing
%             X_rad = [X_rad,X_rad(end)-d_az*[0:NL_stop]]; % Linear
%             X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
%         end
%         
%     else
%         X_rad = [X_rad,d_az*[0:NL_az]]; % Linear
%         X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
%         X_rad = [X_rad,d_az*[NL_az-1:-1:0]]; % Linear
%         X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
%         
%     end
%     X_rad = [X_rad,zeros(1,round(stop_time/prt/2))]; %stop
%     X_rad = [X_rad,zeros(1,round(0.5/prt/2))]; %stop
% end
% X_rad=X_rad(1:Na);


Y_rad=0*(1:epochs);
Z_rad=0*(1:epochs);


r_az = lambda/m_az/max(X_rad);% azimuth resolution!
r_freq = 1/(prt*epochs);% frequency resolution!

%% Targets' Properties:
%%% Position:
% Polar
R_tar = [20,20,20]-2.3; theta_tar = [-10,0,20]; ph_tar=zeros(size(R_tar));
% R_tar = [20]; theta_tar = 1*[0.0]; ph_tar=zeros(size(R_tar));

%%% Displacement functions: 
% Fluctuations:
dR_amp_tar = [0.005,0.001,0.002];%,.00015];
% dR_amp_tar = 1e-3;
dR_f_tar = [7,2,.4];%./prt;%,0.01]; % displacement frequency Hz(1/epochs)
% dR_f_tar = 2;%./prt;%,0.01]; % displacement frequency Hz(1/epochs)


% Cartesian:
X_tar=R_tar.*sind(theta_tar);Y_tar=R_tar.*sqrt(1-sind(theta_tar).^2-sind(ph_tar).^2);Z_tar=R_tar.*sind(ph_tar);
num_tar=length(R_tar);
vib_mode = 'sine'; % Vibration signal model: 'sine' or 'triangle'
if strcmp('SINE',upper(vib_mode))
    dR_tar_fluc = displacement_model_sin_TS( dR_amp_tar,dR_f_tar,prt,time_tot );
elseif strcmp('TRIANGLE',upper(vib_mode))
    dR_tar_fluc = displacement_model_triangle_TS( dR_amp_tar,dR_f_tar,prt,time_tot );
end
% Total displacement:
dR_tar_total_disp = dR_tar_fluc;%+(1e-4)*(rand(Na,length(R_tar))-0.5);%+dR_linear;
figure;
subplot(3,1,1);plot(dR_tar_total_disp,'LineWidth',1.5);xlabel('Time samples');ylabel('Displacement (m)')
subplot(3,1,2);plot(sum(dR_tar_total_disp,2),'LineWidth',1.5);xlabel('Time samples');ylabel('Displacement (m)')
f_dif1 = fft(sum(dR_tar_total_disp,2));

f_ax = r_freq*(0:1:length(f_dif1)-1);
subplot(3,1,3);plot(f_ax,abs(f_dif1));xlabel('f')

if export_flag
    print(gcf,[export_directory , 'Displacements.jpg'],'-djpeg','-r400');
end
%% SIGNAL
TomoSAR_SIGNAL_RT

%% Show parameters in a Table
table_row = {'fc (GHz)';'Wavelength (m)';'Bandwidth (GHz)';'Chirp duration (s)';'snr (dB)';'Range samples (Nr)';'Azimuth samples (Na)';'Azimuth antenna spacing (m)';'Azimuth mode';...
             'Azimuth beam (deg)';'Targets number';'Targets Range (m)';'Targets Phi (deg)';'Targets Theta (deg)';'Targets Vib. Amp (m)';'Targets Vib. f (1/ep)'};
table_val = {fc/1e9;lambda;bw/1e9;T;snr;Nr;Na;d_az;upper(mode_az);beam_az;num_tar;num2str(R_tar);num2str(ph_tar);num2str(theta_tar);num2str(dR_amp_tar);num2str(dR_f_tar)};
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
rc=fft(signal_MIMO_TS,[],1);%rc = rc(1:500,:);
% slc = fft(rc,[181],2);
slc =fftshift(fft(rc,[181],2),2);
[r_tar,c_tar] = find( abs(slc)==max(max(abs(slc))) );
% figure;plot(cumsum(wrapToPi(diff(angle(rc(r_tar,:))))))
figure;imagesc(rad2deg(lambda/2/d_az)*linspace(-.5,.5,181),[],10*log10(abs(slc)));xlabel('Az');ylabel('R');
figure;plot( unwrap(angle(squeeze(rc(r_tar,:))))*lambda/4/pi )
figure;plot( abs(fft(unwrap(angle(squeeze(rc(r_tar,:))))*lambda/4/pi)) )

%% Spectral estimation:
w_r=[r_tar-0,r_tar+0];
antenna_loc = X_rad;
time_vec = (0:epochs-1).*prt;
theta = [-30:2:30]; %round(theta_tar)*[0.5:0.1:1.5]; 
fdisp = [ 0:.2:12 ]; % displacement frequency
Adisp = linspace( 0,8,17)*1e-3; % displacement amplitude
figure;plot(time_vec,antenna_loc,'cyan','LineWidth',2);hold on;scatter(time_vec,antenna_loc,'b','Filled');ylabel('Azimuth locations (m)');xlabel('Time (samples)');title('data acquisition mode')
if export_flag
    print(gcf,[export_directory , 'trajectory.jpg'],'-djpeg','-r400');
end
%
c_iter = 4;
%[Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVib_RT(rc,w_r,lambda,antenna_loc,time_vec,theta,Adisp,fdisp,num_tar,mode_az,vib_mode);
[Pbf,Pcapon,Pmusic,Pbf_C,Pcapon_C,Pmusic_C]=SpectralEstimation_AzVib_RT_CLEAN(rc,w_r,lambda,antenna_loc,time_vec,theta,Adisp,fdisp,num_tar,mode_az,vib_mode,c_iter);
%
Pbf = norm01((Pbf)); Pcapon = norm01((Pcapon)); Pmusic = norm01((Pmusic)); %Pmf = norm01((Pmf)); 
Pbf_C = norm01((Pbf_C)); Pcapon_C = norm01((Pcapon_C)); Pmusic_C = norm01((Pmusic_C)); %Pmf = norm01((Pmf)); 

Pbf = absdB(Pbf); Pcapon = absdB(Pcapon); Pmusic = absdB(Pmusic);
Pbf_C = absdB(Pbf_C); Pcapon_C = absdB(Pcapon_C); Pmusic_C = absdB(Pmusic_C);
caxislim=[-10,0]; 

%% RESULTS:
[rr,cc,vv]=ind2sub(size(Pbf),find(Pbf==max(Pbf,[],'all')));          thet_bf    = theta(rr),Ampd_bf   =Adisp(cc), freq_bf    = fdisp(vv)
[rr,cc,vv]=ind2sub(size(Pcapon),find(Pcapon==max(Pcapon,[],'all'))); thet_capon = theta(rr),Ampd_capon=Adisp(cc), freq_capon = fdisp(vv)
[rr,cc,vv]=ind2sub(size(Pmusic),find(Pmusic==max(Pmusic,[],'all'))); thet_music = theta(rr),Ampd_music=Adisp(cc), freq_music = fdisp(vv)

snr_bf    = 10*log10( 1/std(Pbf(15:end,15:end,8:end),[],'all') );
snr_capon = 10*log10( 1/std(Pcapon(15:end,15:end,8:end),[],'all') );
snr_music = 10*log10( 1/std(Pmusic(15:end,15:end,8:end),[],'all') );

for tar_qq = 1:num_tar
    thet_idx (tar_qq) = find(theta == theta_tar (tar_qq));
    Adisp_idx(tar_qq) = find(Adisp == dR_amp_tar(tar_qq));
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

