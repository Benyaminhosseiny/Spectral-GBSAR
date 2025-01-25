% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
function [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVel_RT(rc,w_r,lambda,antenna_loc,time_vec,theta,vel,num_tar,mode)
% rc: range compressed data-cube: range*azimuth*elevation
% w_r: neighborhood samples in range direction for estimating covariance matrix[2*1 vector array]. 
% lambda: wavelength (m)
% antenna_loc: SAR antenna locations (1d array)
% time_vec: Data acquisition times (1d array)
% theta (deg): AOA search span (1d array)
% Adisp (m): fluctuation amplitude search span (1d array)
% fdisp (1/s): fluctuation frequency search span (1d array)
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
time_vec = time_vec(:);
antenna_loc = antenna_loc(:);

%% 1-Covariance
x = squeeze(rc(w_r(1):w_r(2),:)); %x =reshape(x,[],size(x,3));
Rx0 = cov(x);
if length(Rx0)==1
    Rx0 = x'*x;
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

URA_az_vel_RT = @(theta,vel)exp( -1j*2*pi*( m*antenna_loc*sind(theta) + 2*time_vec*vel )/lambda );  % Steering vector: + displacement frequency (fluctuations)

for theta_ii = 1:length(theta)
    for vel_ii = 1:length(vel)
            SS = URA_az_vel_RT(theta(theta_ii), vel(vel_ii)); % Fluctuations
            
            % BF:
            PP = SS'*Rx*SS; % pseudospectrum estimate
            Pbf(theta_ii,vel_ii) = PP;
            
            % Capon:
%             PP = SS'*Rx_inv*SS; % pseudospectrum estimate
            PP = SS'/Rx*SS; % pseudospectrum estimate
            Pcapon(theta_ii,vel_ii) = 1/PP;
            
            % MUSIC:
            PP = SS'*(Vn*Vn')*SS; % pseudospectrum estimate
            Pmusic(theta_ii,vel_ii) = 1/ PP;
    end
end

Pbf = abs(Pbf); 
Pcapon = abs(Pcapon);
Pmusic = abs(Pmusic); 

end