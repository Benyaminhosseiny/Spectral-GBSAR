% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
clear;clc;close all
addpath('./SES_functions')
c = physconst('LightSpeed');
export_flag = 0;
export_directory = '';
if export_flag
    mkdir(export_directory)
    export_directory=[export_directory,'Sim-AzVel\',datestr(now,'yyyymmdd_HHMM'),'\'];
    mkdir(export_directory);
end
fontsizefig = 14; fontname = 'times'; % Set it to times
set(0,'DefaultAxesFontName',fontname,'DefaultTextFontName',fontname,'defaultAxesFontSize',fontsizefig);

%% Parameters:
fc = 77e9; lambda = c/fc;
bw = 1e9;
T = 60e-6; cr = bw/T;
snr = -5;
Nr = 512;
Na = 8;
Ne = 1;
epochs = 50;
prt = 0.02; % (s)
d_az = lambda/2;
d_el = lambda/2;
time_tot = prt*(epochs-1); % Total monitoring time

signal_model = 'polar'; % 'cartesian' or 'polar'
mode_az='mimo';
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
X_rad=d_az*(0:Na-1);
Y_rad=0;
Z_rad=0*d_az*(0:epochs-1);

%% Target position:
% Polar:
R_tar = [20,20,20]; theta_tar = [-20,10,0]; ph_tar=zeros(size(R_tar));

% Cartesian:
X_tar=R_tar.*sind(theta_tar);Y_tar=R_tar.*sqrt(1-sind(theta_tar).^2-sind(ph_tar).^2);Z_tar=R_tar.*sind(ph_tar);
num_tar=length(R_tar);

%% Displacement functions: 
% Linear displacements:
vel_tar = [-0.015,0.005,0]; % m/s
dR_linear = (0:prt:time_tot)'*vel_tar;

% Total displacement:
dR_tar_total_disp = dR_linear;
figure;plot(dR_tar_total_disp,'LineWidth',1.5);xlabel('Time samples');ylabel('Displacement (m)')

%% SIGNAL:
TomoSAR_SIGNAL

%% Show parameters in a Table
table_row = {'fc (GHz)';'Wavelength (m)';'Bandwidth (GHz)';'Chirp duration (s)';'snr (dB)';'Range samples (Nr)';'Azimuth samples (Na)';'Azimuth antenna spacing (m)';...
             'Azimuth mode';'Azimuth beam (deg)';'epochs';'Targets number';'Targets Range (m)';'Targets Phi (deg)';'Targets Theta (deg)';'Targets velocity (m/ep)'};
table_val = {fc/1e9;lambda;bw/1e9;T;snr;Nr;Na;d_az;upper(mode_az);beam_az;epochs;num_tar;num2str(R_tar);num2str(ph_tar);num2str(theta_tar);num2str(vel_tar)};
Table = table(table_val,'RowNames',table_row);
figure;tbl=uitable('Data', table_val, 'RowName', table_row,'ColumnName',{'Value'},'Position',[10,10,500,400]);tbl.ColumnWidth={150};
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
figure;subplot(1,2,1);imagesc(abs(slc(:,:,1)));xlabel('Az');ylabel('R');

%% Spectral estimation:
w_r   = [r_tar-0,r_tar+0];
theta = (-26:2:26);-round(beam_az/2):1:round(beam_az/2); 
R0 = r_tar*c/2/bw;
ant_loc_az = d_az*(0:Na-1)';
time_vec = (0:epochs-1).*prt;
% vel_res = lambda/2/(epochs-1);
vel=(-0.03:0.001:0.03);

% Random sampling for Non-Uniform baselines:
% rs = sort(randperm(epochs,round(epochs/4)));time_vec=time_vec(rs);time_vec=time_vec-min(time_vec);rc=rc(:,:,rs);

%
% [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVel(rc,w_r,lambda,ant_loc_az,time_vec,theta,vel,1,mode_az);
c_iter = 4;
[Pbf,Pcapon,Pmusic,Pbf_C,Pcapon_C,Pmusic_C]=SpectralEstimation_AzVel_CLEAN(rc,w_r,lambda,ant_loc_az,time_vec,theta,vel,num_tar,mode_az,c_iter);

%
slc_vel =fftshift(fft(slc,length(vel),3),3);
vel_ax_full = linspace(-lambda/4,lambda/4,size(slc_vel,3))./prt;
theta_ax_full=linspace(-beam_az/2,beam_az/2,size(slc_vel,2)); 
Pft = abs((squeeze(slc_vel(r_tar,:,:))));
Pft = norm01((Pft)); Pbf = norm01((Pbf)); Pcapon = norm01((Pcapon)); Pmusic = norm01((Pmusic));
Pbf_C = norm01((Pbf_C)); Pcapon_C = norm01((Pcapon_C)); Pmusic_C = norm01((Pmusic_C));
Pft = absdB(Pft);    Pbf = absdB(Pbf);    Pcapon = absdB(Pcapon);    Pmusic = absdB(Pmusic);
Pbf_C = absdB(Pbf_C); Pcapon_C = absdB(Pcapon_C);    Pmusic_C = absdB(Pmusic_C);
caxislim=[-10,0]; 

% % %% RESULTS:
[rrr,cc]=find(Pbf==1);   veloc_bf    = vel(cc); ph_bf   =theta(rrr);
[rrr,cc]=find(Pcapon==1);veloc_capon = vel(cc); ph_capon=theta(rrr);
[rrr,cc]=find(Pmusic==1);veloc_music = vel(cc); ph_music=theta(rrr);



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
subplot(1,4,4); imagesc(vel_ax_full,theta_ax_full,Pft);xlim([min(vel),max(vel)]);ylim([min(theta),max(theta)]);%cbar.Title.String = "[dB]";
hold on; scatter(vel_tar',theta_tar,200,'ro','LineWidth',1.5);xlabel([{'{\it v} (m/s)'}]);
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
