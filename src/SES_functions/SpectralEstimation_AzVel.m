% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
function [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVel(rc,w_r,lambda,ant_loc_az,time_vec,theta,vel,num_tar,mode)
% rc: range compressed data-cube: range*azimuth*elevation
% w_r: neighborhood samples in range direction for estimating covariance matrix[2*1 vector array]. 
% lambda: wavelength (m)
% ant_loc_az: SAR antenna locations in Azimuth direction (1d array)
% time_vec: Data acquisition times (1d array)
% theta (deg): AOA search span in Azimuth direction (1d array)
% vel (m/s): linear velocity search span (1d array)
% num_tar: Number of targets (for MUSIC)
% mode: SAR imaging mode: 'mono' or 'mimo' (default: 'mimo')
if nargin<9
    mode ='mimo';
end

if strcmp('MIMO',upper(mode))
    m=1;
elseif strcmp('MONO',upper(mode))
    m=2;
end
%
ant_loc_az = ant_loc_az(:);
time_vec = time_vec(:)';

%% 1-Covariance
x = rc(w_r(1):w_r(2),:,:); x =reshape(x,size(x,1),[]);
if size(x,1)>1
    Rx0 = cov(x);
else
    Rx0=x'*x;
end
Rx = Rx0;
% Normalize Covariance Matrix, where the diagonal items are 1 (better results on BF and Capone)
E = diag(diag(Rx));
Rx = E^(-1/2)*Rx*E^(-1/2);
Rx = Rx + 1e-3*eye(length(Rx));
% figure; sgtitle('Cov');subplot(1,2,1);imagesc(abs(Rx));title('abs');subplot(1,2,2);imagesc(angle(Rx));title('angle');

%% 2-Spectral estimation:
Rx_inv = inv(Rx);
[eigenVec,eigenVal] = eig(Rx0); 
Vn = eigenVec(:,1:end-num_tar); % Noise Subspace

URA_Az_Vel = @(theta,vel)exp( -1j*2*pi*( m*ant_loc_az*sind(theta)+2*time_vec*vel )/lambda );  % Steering vector Az-Velocity

Pbf=zeros([length(theta),length(vel)]);
Pcapon=zeros([length(theta),length(vel)]);
Pmusic=zeros([length(theta),length(vel)]);
for thet_ii = 1:length(theta)
    for vel_ii = 1:length(vel)
        % Steering (Sensing) vector:
        SS = URA_Az_Vel(theta(thet_ii),vel(vel_ii));
        SS = SS(:);
        
        % BF:
        PP = SS'*Rx*SS;
        Pbf(thet_ii,vel_ii) = PP; 
        %Capon:
%         PP = SS'*Rx_inv*SS;
        PP = SS'/Rx*SS;
        Pcapon(thet_ii,vel_ii) = 1/PP; 
        % MUSIC:
        PP = SS'*(Vn*Vn')*SS; 
        Pmusic(thet_ii,vel_ii) = 1/ PP; 
    end
end

Pbf = abs(Pbf); 
Pcapon = abs(Pcapon);
Pmusic = abs(Pmusic); 

end