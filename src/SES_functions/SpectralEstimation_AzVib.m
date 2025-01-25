% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
function [Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVib(rc,w_r,lambda,time_vec,antenna_loc,theta,Adisp,fdisp,num_tar,mode,vib_mode)
% rc: range compressed data-cube: range*azimuth*elevation
% w_r: neighborhood samples in range direction for estimating covariance matrix[2*1 vector array]. 
% lambda: wavelength (m)
% time_vec: Data acquisition times (1d array)
% Adisp (m): fluctuation amplitude search span (1d array)
% fdisp (1/s): fluctuation frequency search span (1d array)
% num_tar: Number of targets (for MUSIC)
% mode: SAR imaging mode: 'mono' or 'mimo' (default: 'mimo')
if nargin<10
    mode ='mimo';
end
if nargin<11
    vib_mode='SINE';
end
if strcmp('MIMO',upper(mode))
    m=1;
elseif strcmp('MONO',upper(mode))
    m=2;
end
%
antenna_loc = antenna_loc(:);
time_vec = time_vec(:)';

% ph0 = linspace(0,2*pi,10);
ph0 = 0;

%% 1-Covariance
x = rc(w_r(1):w_r(2),:,:); x =reshape(x,size(x,1),[]);
% Rx0 = cov(x);
Rx0 = x'*x;
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

if strcmp('SINE',upper(vib_mode))
    URA_az_vib_fluc = @(theta,Adisp,fdisp,ph0)exp( -1j*2*pi*(m*antenna_loc*sind(theta) + 2*Adisp*sin(2*pi*fdisp*time_vec+ph0))/lambda );  % Steering vector: linear velocity + displacement frequency (fluctuations)
elseif strcmp('TRIANGLE',upper(vib_mode))
    Nh = 3;7;5;         % number of harmonics [the higher the Nh the better approximation of tiangular signal]
    ii = (0:Nh-1)'; % harmonic label
    nn = 2*ii+1;    % harmonic mode number
    URA_az_vib_fluc = @(theta,Adisp,fdisp,ph0)exp( -1j*2*pi*(m*antenna_loc*sind(theta) + 2*Adisp*sum(8/(pi^2)*((-1).^ii).*(nn.^-2).*sin(2*pi*fdisp*nn*time_vec + ph0),1))/lambda );  % Steering vector: linear velocity + displacement frequency (fluctuations)
end

Pbf    = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));
Pcapon = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));
Pmusic = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));
for theta_ii = 1:length(theta)
    for fd_ii = 1:length(fdisp)
        for Ad_ii=1:length(Adisp)
            for ph0_ii = 1:length(ph0)
                SS = URA_az_vib_fluc(theta(theta_ii), Adisp(Ad_ii), fdisp(fd_ii),ph0(ph0_ii)); % Fluctuations
                SS = SS(:);
                
                % BF:
                PP = SS'*Rx0*SS; % pseudospectrum estimate
                Pbf(theta_ii,Ad_ii,fd_ii,ph0_ii) = PP;
%                 
%                 %Capon:
                PP = SS'*Rx_inv*SS; % pseudospectrum estimate
                Pcapon(theta_ii,Ad_ii,fd_ii,ph0_ii) = 1/PP;
%                 
                %MUSIC
                PP = SS'*(Vn*Vn')*SS; % pseudospectrum estimate
                Pmusic(theta_ii,Ad_ii,fd_ii,ph0_ii) = 1/ PP;
            end
        end
    end
end

Pbf = abs(Pbf);
Pcapon = abs(Pcapon);
Pmusic = abs(Pmusic);


Pbf=squeeze(max(Pbf,[],4));
Pcapon=squeeze(max(Pcapon,[],4));
Pmusic=squeeze(max(Pmusic,[],4));

end