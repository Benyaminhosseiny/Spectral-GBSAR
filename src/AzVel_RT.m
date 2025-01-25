% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
clear;clc;
% close all
addpath('./SES_functions')
export_flag = 0;
export_directory = '';
if export_flag
    mkdir(export_directory)
    export_directory=[export_directory,'Sim-AzVel_RT\',datestr(now,'yyyymmdd_HHMM'),'\'];
    mkdir(export_directory);
end
fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);

%% Parameters:
c = physconst('LightSpeed');
fc = 77e9; lambda = c/fc;
bw = 1e9/2;
T = 60e-6; cr = bw/T;
snr = 5;
Nr = 512;
Na = 300;
Ne = 1;
epochs = Na;
prt = 0.02; % (s)
radar_vel=0.025;0.0245; %m/s
d_az = prt*radar_vel;%4*lambda/16; %antenna spacing
d_az = lambda/8;
d_el = lambda/4;
time_tot = prt*(epochs-1); % Total monitoring time

signal_model = 'polar'; % 'cartesian' or 'polar'
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
% X_rad=d_az*(0:Na-1);
% Random scenario:
% L_az = 150e-3;
% X_rad=linspace(0,L_az,Na);X_rad = 1*[ 0,X_rad(randperm(length(X_rad)-1)) ];

% Zigzag scenario:
L_az = 1.*40e-3; NL_az = round(L_az/d_az);
L_stop = 1*40e-3; NL_stop = round(L_stop/d_az);
stop_time=0*0.5; % s
X_rad=[0];
for ii=1:floor(epochs/NL_az/2)+1
    if stop_time*L_stop~=0
        for jj=1:floor(L_az/L_stop) % increasing
            X_rad = [X_rad,X_rad(end)+d_az*[0:NL_stop]]; % Linear
            X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
        end
        X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt/2))]; %stop
        for jj=1:floor(L_az/L_stop) % decreasing
            X_rad = [X_rad,X_rad(end)-d_az*[0:NL_stop]]; % Linear
            X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
        end
    else
        X_rad = [X_rad,d_az*[0:NL_az]]; % Linear
        X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
        X_rad = [X_rad,d_az*[NL_az-1:-1:0]]; % Linear
        X_rad = [X_rad,X_rad(end)*ones(1,round(stop_time/prt))]; %stop
        
    end
    X_rad = [X_rad,zeros(1,round(stop_time/prt/2))]; %stop
%     X_rad = [X_rad,zeros(1,round(0.5/prt/2))]; %stop
end
X_rad=X_rad(1:Na);

Y_rad=0*(1:epochs);
Z_rad=0*(1:epochs);

%% Target position:
% Polar
R_tar = [20,20,20]; theta_tar = [-20,10,0]; ph_tar=zeros(size(R_tar));

% Cartesian:
X_tar=R_tar.*sind(theta_tar);Y_tar=R_tar.*sqrt(1-sind(theta_tar).^2-sind(ph_tar).^2);Z_tar=R_tar.*sind(ph_tar);
num_tar=length(R_tar);

%% Displacement functions:
% Linear displacement:
vel_tar = [-0.0006,0.0002,0]./prt; % m/s
vel_tar = [-0.015,0.005,0]; % m/s
dR_linear = (0:prt:time_tot)'*vel_tar;
% Total displacement:
dR_tar_total_disp = dR_linear;
figure;plot(dR_tar_total_disp,'LineWidth',1.5);xlabel('Time samples');ylabel('Displacement (m)')

%% SIGNAL
SIGNAL_RT

%% Show parameters in a Table
table_row = {'fc (GHz)';'Wavelength (m)';'Bandwidth (GHz)';'Chirp duration (s)';'snr (dB)';'Range samples (Nr)';'Azimuth samples (Na)';'SAR length (m)';'Azimuth antenna spacing (m)';...
    'Azimuth mode';'Azimuth beam (deg)';'Targets number';'Targets Range (m)';'Targets Phi (deg)';'Targets Theta (deg)';'Targets velocity (m/ep)'};
table_val = {fc/1e9;lambda;bw/1e9;T;snr;Nr;Na;L_az;d_az;upper(mode_az);beam_az;num_tar;num2str(R_tar);num2str(ph_tar);num2str(theta_tar);num2str(vel_tar)};
Table = table(table_val,'RowNames',table_row);
figure;tbl=uitable('Data', table_val, 'RowName', table_row,'ColumnName',{'Value'},'Position',[10,10,500,400]);tbl.ColumnWidth={10+num_tar*50};
if export_flag
    print(gcf,[export_directory , 'Table of Parameters.jpg'],'-djpeg','-r400');
end
%% ========================================================================
%% +-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*
%% +-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*+-*
%% ========================================================================

%% Range compression:
rc=fft(signal_MIMO_TS,[],1);
slc =fftshift(fft(rc,[1*Na],2),2);
[r_tar,c_tar] = find( abs(slc)==max(max(abs(slc))) );
figure;imagesc(abs(slc));xlabel('Az');ylabel('R');
figure;plot( unwrap(angle(squeeze(rc(r_tar,:))))*lambda/4/pi )



%% Spectral estimation:
w_r   = [r_tar-2,r_tar+2];
w_r   = [r_tar,r_tar];
antenna_loc = X_rad;
time_vec = (0:Na-1).*prt;
theta = -26:2:26;%-round(beam_az/16):2:round(beam_az/16);
vel=(-0.03:0.001:0.03);

% Random sampling for Non-Uniform baselines:
% rs = sort(randperm(epochs,round(epochs/5)));time_vec=time_vec(rs);antenna_loc=antenna_loc(rs);rc=rc(:,rs);

figure('Position', [30 60 1000 500]);plot(time_vec,antenna_loc,'cyan','LineWidth',2);hold on;scatter(time_vec,antenna_loc,10,'b','Filled');ylabel('Azimuth locations (m)');xlabel('Time (samples)');title('data acquisition mode')
if export_flag
    print(gcf,[export_directory , 'Baseline distribution.jpg'],'-djpeg','-r400');
end
c_iter = 4;
% [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVel_RT(rc,w_r,lambda,antenna_loc,time_vec,theta,vel,num_tar,mode_az);
[Pbf, Pcapon, Pmusic, Pbf_C, Pcapon_C, Pmusic_C] = SpectralEstimation_AzVel_RT_CLEAN(rc,w_r,lambda,antenna_loc,time_vec,theta,vel,num_tar,mode_az,c_iter);
%
Pbf = norm01(Pbf); Pcapon = norm01(Pcapon); Pmusic = norm01(Pmusic);
Pbf_C = norm01(Pbf_C); Pcapon_C = norm01(Pcapon_C); Pmusic_C = norm01(Pmusic_C);
Pbf = absdB(Pbf); Pcapon = absdB(Pcapon); Pmusic = absdB(Pmusic);
Pbf_C = absdB(Pbf_C); Pcapon_C = absdB(Pcapon_C); Pmusic_C = absdB(Pmusic_C);
caxislim=[-10,0]; 
%% RESULTS
[rr,cc]=find(Pbf==max(Pbf,[],'all'));      veloc_bf    = vel(cc), thet_bf   =theta(rr)
[rr,cc]=find(Pcapon==max(Pcapon,[],'all'));veloc_capon = vel(cc), thet_capon=theta(rr)
[rr,cc]=find(Pmusic==max(Pmusic,[],'all'));veloc_music = vel(cc), thet_music=theta(rr)

for tar_qq = 1:num_tar
    thet_idx (tar_qq) = find(theta == theta_tar (tar_qq));
    vel_idx(tar_qq) = find(round(vel,3) == round(vel_tar  (tar_qq),3));
end
P_idx = sub2ind([length(theta), length(vel)],thet_idx,vel_idx);
snr_bf      = 10*log10( mean(Pbf     (P_idx))/((sum( Pbf     (:) )-sum( Pbf     ( P_idx ) ))/(length(Pbf     (:))-num_tar)) );
snr_capon   = 10*log10( mean(Pcapon  (P_idx))/((sum( Pcapon  (:) )-sum( Pcapon  ( P_idx ) ))/(length(Pcapon  (:))-num_tar)) );
snr_music   = 10*log10( mean(Pmusic  (P_idx))/((sum( Pmusic  (:) )-sum( Pmusic  ( P_idx ) ))/(length(Pmusic  (:))-num_tar)) );
snr_bf_C    = 10*log10( mean(Pbf_C   (P_idx))/((sum( Pbf_C   (:) )-sum( Pbf_C   ( P_idx ) ))/(length(Pbf_C   (:))-num_tar)) );
snr_capon_C = 10*log10( mean(Pcapon_C(P_idx))/((sum( Pcapon_C(:) )-sum( Pcapon_C( P_idx ) ))/(length(Pcapon_C(:))-num_tar)) );
snr_music_C = 10*log10( mean(Pmusic_C(P_idx))/((sum( Pmusic_C(:) )-sum( Pmusic_C( P_idx ) ))/(length(Pmusic_C(:))-num_tar)) );

set(0,'DefaultAxesFontName','times','DefaultTextFontName','times','defaultAxesFontSize',12);
figure('Position', [30 60 2000 400]); sgtitle('Azimuth/Velocity')
subplot(1,4,1); imagesc(vel,theta,Pbf);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_bf)]]);
subplot(1,4,2); imagesc(vel,theta,Pcapon);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_capon)]]);
subplot(1,4,3); imagesc(vel,theta,Pmusic);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_music)]]);
cbar = colorbar('Position',[0.93 0.168 0.022 0.7]);caxis(caxislim)
cbar.Title.String = "[dB]";
Lgnd = legend({'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.92; Lgnd.Position(2) = 0.01;
if export_flag
    print(gcf,[export_directory , '2D normalized power spectrum density.jpg'],'-djpeg','-r400');
end
% CLEAN:
set(0,'DefaultAxesFontName','times','DefaultTextFontName','times','defaultAxesFontSize',12);
figure('Position', [30 60 2000 400]); sgtitle('Azimuth/Velocity')
subplot(1,4,1); imagesc(vel,theta,Pbf_C);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_bf_C)]]);
subplot(1,4,2); imagesc(vel,theta,Pcapon_C);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_capon_C)]]);
subplot(1,4,3); imagesc(vel,theta,Pmusic_C);caxis(caxislim)
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);ylabel('\theta (deg)');xlabel([{'{\it v} (m/s)'},['SNR ', num2str(snr_music_C)]]);
cbar = colorbar('Position',[0.93 0.168 0.022 0.7]);caxis(caxislim);%caxis([0,1])
cbar.Title.String = "[dB]";
Lgnd = legend({'Ground Truth'},'FontSize',10); Lgnd.Position(1) = 0.92; Lgnd.Position(2) = 0.01;
if export_flag
    print(gcf,[export_directory , '2D normalized power spectrum density_CLEAN.jpg'],'-djpeg','-r400');
end


for ii = 1:num_tar
    [~,theta_tar_idx(ii)]=min(abs(theta-theta_tar(ii)));
    [~,vel_tar_idx(ii)]=min(abs(vel-vel_tar(ii)));
end
w_snr=[3,7];% 3+1+3=7 =>7*7 w
for ii = 1:num_tar
    idx_r = [max([1, theta_tar_idx(ii)-w_snr(1)]), min([length(theta), theta_tar_idx(ii)+w_snr(2)])];
    idx_c = [max([1, vel_tar_idx(ii)-w_snr(2)]), min([length(vel), vel_tar_idx(ii)+w_snr(2)])];
    P_tar_bf(ii) = Pbf(theta_tar_idx(ii),vel_tar_idx(ii));
    P_tar_capon(ii) = Pcapon(theta_tar_idx(ii),vel_tar_idx(ii));
    P_tar_music(ii) = Pmusic(theta_tar_idx(ii),vel_tar_idx(ii));
    
    snr_bf(ii)    = 10*log10( SNR_Custom(Pbf(idx_r(1):idx_r(2),idx_c(1):idx_c(2))) );
    snr_capon(ii) = 10*log10( SNR_Custom(Pcapon(idx_r(1):idx_r(2),idx_c(1):idx_c(2))) );
    snr_music(ii) = 10*log10( SNR_Custom(Pmusic(idx_r(1):idx_r(2),idx_c(1):idx_c(2))) );
    
    pslr_bf(ii)    = PSLR( Pbf(idx_r(1):idx_r(2),idx_c(1):idx_c(2)) )/2;
    pslr_capon(ii) = PSLR( Pcapon(idx_r(1):idx_r(2),idx_c(1):idx_c(2)) )/2;
    pslr_music(ii) = PSLR( Pmusic(idx_r(1):idx_r(2),idx_c(1):idx_c(2)) )/2;
end

table_row = {'Aperture length (m)'; 'Aperture step length (m)'; 'Stop time (s)'; 'Theta (deg)'; 'v(m/s)'; 'BF';'PSD BF';'SNR (dB) BF';'PSLR (dB) BF';'CAPON';'PSD CAPON';'SNR (dB) CAPON';'PSLR (dB) CAPON';'MUSIC';'PSD MUSIC';'SNR (dB) MUSIC';'PSLR (dB) MUSIC'};
table_val = {num2str(L_az); num2str(L_stop); num2str(stop_time); num2str(theta_tar); num2str(vel_tar); ''; num2str(P_tar_bf); num2str(snr_bf); num2str(pslr_bf); ''; num2str(P_tar_capon); ...
    num2str(snr_capon); num2str(pslr_capon); ''; num2str(P_tar_music); num2str(snr_music); num2str(pslr_music)};
Table = table(table_val,'RowNames',table_row);
figure;tbl=uitable('Data', table_val, 'RowName', table_row,'ColumnName',{'Value'},'Position',[10,10,500,400]);tbl.ColumnWidth={10+num_tar*50};
if export_flag
    print(gcf,[export_directory , 'Table of Evaluations.jpg'],'-djpeg','-r400');
end