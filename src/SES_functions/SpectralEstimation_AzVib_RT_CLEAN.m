% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
function [Pbf0,Pcapon0,Pmusic0,PbfC,PcaponC,PmusicC]=SpectralEstimation_AzVib_RT_CLEAN(rc,w_r,lambda,antenna_loc,time_vec,theta,Adisp,fdisp,num_tar,mode,vib_mode,c_iter)
% * * * By Benyamin Hosseiny * * *
%   * Latest Update: 2023-09-12 *

% rc: range compressed data-cube: range*azimuth*elevation
% w_r: neighborhood samples in range direction for estimating covariance matrix[2*1 vector array].
% Rx: Covariance matrix
% lambda: wavelength (m)
% antenna_loc: SAR antenna locations (1d array)
% time_vec: Data acquisition times (1d array)
% theta (deg): AOA search span (1d array)
% Adisp (m): fluctuation amplitude search span (1d array)
% fdisp (1/s): fluctuation frequency search span (1d array)
% num_tar: Number of targets (for MUSIC)
% mode: SAR imaging mode: 'mono' or 'mimo' (default: 'mimo')
% vib_mode: Vibration model: 'SINE' or 'TRIANGLE' (default: 'SINE')
% c_iter: number of iterations for CLEAN algorithm

if nargin<10
    mode ='mimo';
end
if nargin<11
    vib_mode='SINE';
end
if nargin<12
    c_iter=1;
end

if strcmp('MIMO',upper(mode))
    m=1;
elseif strcmp('MONO',upper(mode))
    m=2;
end
%
time_vec = time_vec(:);
antenna_loc = antenna_loc(:);

%% 1-Covariance:
x = squeeze(rc(w_r(1):w_r(2),:));
if size(x,1)>1
    Rx0 = cov(x);
else
    Rx0 = x'*x;
end
%
Rx = Rx0;
% Normalize Covariance Matrix, where the diagonal items are 1 (better results on BF and Capone)
E = diag(diag(Rx));
Rx = E^(-1/2)*Rx*E^(-1/2);
Rx = Rx + 1e-3*eye(length(Rx));
% figure; sgtitle('Cov');subplot(1,2,1);imagesc(abs(Rx));title('abs');subplot(1,2,2);imagesc(angle(Rx));title('angle');

%% 2-Spectral estimation:
% Rx_inv = inv(Rx);

% ph0 = linspace(0,2*pi,10);
ph0 = 0; 


% Steering vector: + displacement frequency (fluctuations):
if strcmp('SINE',upper(vib_mode))
    URA_az_fluc_RT = @(theta,Adisp,fdisp,ph0)exp( -1j*2*pi*( m*antenna_loc*sind(theta) + 2*Adisp*sin(2*pi*fdisp*time_vec + ph0) )/lambda ); % MIMO
elseif strcmp('TRIANGLE',upper(vib_mode))
    Nh = 2;1;3;7;4;5;         % number of harmonics [the higher the Nh the better approximation of tiangular signal]
    ii = (0:Nh-1)'; % harmonic label
    nn = 2*ii+1;    % harmonic mode number
    URA_az_fluc_RT = @(theta,Adisp,fdisp,ph0)exp( -1j*2*pi*( m*antenna_loc*sind(theta) + 2*Adisp*sum(8/(pi^2)*((-1).^ii).*(nn.^-2).*sin(2*pi*fdisp*nn*time_vec' + ph0),1)')/lambda );  % Steering vector: linear velocity + displacement frequency (fluctuations)
end

PbfC    = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));
PcaponC = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));
PmusicC = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));

Pbf_clean    = zeros(length(theta),length(Adisp),length(fdisp),length(ph0)); 
Pmusic_clean = zeros(length(theta),length(Adisp),length(fdisp),length(ph0)); 
Pcapon_clean = zeros(length(theta),length(Adisp),length(fdisp),length(ph0));

%Initialise trimmed cross spectral matrix (CSM) by setting the diagonal to zero
Rx0(logical(eye(size(Rx0,1)))) = 0;
Rx(logical(eye(size(Rx,1)))) = 0;

%Initialise break criterion
sumOfCSM = sum(sum(abs(Rx0)));
sumOfDegradedCSM = sumOfCSM;

Rx_M = Rx0;
% Rx_M = Rx;
Rxbf = Rx0;
% Rxbf = Rx;
[eigenVec,eigenVal] = eig(Rx_M);
Vn = eigenVec(:,1:end-1); % Noise Subspace
    
maxIterations = c_iter; loopGain = 0.9; % CLEAN parameters
for cleanMapIterations = 1:maxIterations
    Rx_inv = inv(Rx);
    if cleanMapIterations>1
        [eigenVec,eigenVal] = eig(Rx_M);
        Vn = eigenVec(:,1:end-1);
    end

    for theta_ii = 1:length(theta)
        for Ad_ii=1:length(Adisp)
            for fd_ii = 1:length(fdisp)
                for ph0_ii = 1:length(ph0)
                    SS = URA_az_fluc_RT( theta(theta_ii), Adisp(Ad_ii), fdisp(fd_ii), ph0(ph0_ii) ); % Fluctuations

                    % BF:
                    PP = SS'*Rxbf*SS; % pseudospectrum estimate
                    PbfC(theta_ii,Ad_ii,fd_ii,ph0_ii) = PP;

                    % Capon:
                    PP = SS'*Rx_inv*SS; % pseudospectrum estimate
                    % PP = SS'/Rx*SS; % pseudospectrum estimate
                    PcaponC(theta_ii,Ad_ii,fd_ii,ph0_ii) = 1/PP;
                    
                    
                    % MUSIC:
                    PP = SS'*(Vn*Vn')*SS; % pseudospectrum estimate
%                     PP = SS'*Rx_M*SS; % pseudospectrum estimate
                    PmusicC(theta_ii,Ad_ii,fd_ii,ph0_ii) = 1/ PP;

                end
            end
        end
    end
    if cleanMapIterations == 1
        Pbf0    = abs(PbfC);
        Pcapon0 = abs(PcaponC);
        Pmusic0 = abs(PmusicC);
    end
    %% CLEAN:
    % -------------------------------------------------------
    % 2. Find peak value and its position in dirty map
%     [maxPeakValue, maxPeakIndx] = max( abs( Pbf(:) ) );
%     [maxPeakValueIndx_thet, maxPeakValueIndx_A, maxPeakValueIndx_f] = ind2sub( size(Pbf), maxPeakIndx );


    % Find peaks BF
    [peaks, locations] = findpeaks( abs( PbfC(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pbf, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pbf = locations(sortedIndices);
    
    % Find peaks Capon
    [peaks, locations] = findpeaks( abs( PcaponC(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pcapon, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pcapon = locations(sortedIndices);
    
    % Find peaks MUSIC
    [peaks, locations] = findpeaks( abs( PmusicC(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pmusic, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pmusic = locations(sortedIndices);

    % Clean beam with specified width and max value of 1
    PmaxCleanBeam = zeros(length(theta), length(Adisp), length(fdisp));
    g = 0;
    for tt = 1:num_tar % Three scatterers
        [maxPeakValueIndx_thet(tt), maxPeakValueIndx_A(tt), maxPeakValueIndx_f(tt)] = ind2sub( size(PbfC), sortedLocations_Pbf(tt) );
        
        % -------------------------------------------------------
        % 3. Calculate the CSM induced by the peak source
        
        % Steering vector to location of peak source
        g = 0+URA_az_fluc_RT(theta(maxPeakValueIndx_thet(tt)), Adisp(maxPeakValueIndx_A(tt)), fdisp(maxPeakValueIndx_f(tt)), ph0(ph0_ii));
                
        % Cross spectral matrix induced by peak source in that direction (eq. 11)
        G = g*g';
        G(logical(eye(size(G,1)))) = 0;
        
        
        
        % -------------------------------------------------------
        % 4. New updated map with clean beam from peak source location
        % Clean beam with specified width and max value of 1
        
        PmaxCleanBeam(maxPeakValueIndx_thet(tt), maxPeakValueIndx_A(tt), maxPeakValueIndx_f(tt)) = 1;
        
        % Update clean map with clean beam from peak source location
        Pbf_clean    = Pbf_clean    + loopGain*maxPeakValue_Pbf(tt)   *PmaxCleanBeam;
        Pcapon_clean = Pcapon_clean + loopGain*maxPeakValue_Pcapon(tt)*PmaxCleanBeam;
        Pmusic_clean = Pmusic_clean + loopGain*maxPeakValue_Pmusic(tt)*PmaxCleanBeam;
        
        % -------------------------------------------------------
        % 5. Calculate degraded cross spectral matrix
        % Basically removing the PSF from that location of the plot
        Rxbf = Rxbf - loopGain*maxPeakValue_Pbf(tt)*G;
        Rx   = Rx   - loopGain*maxPeakValue_Pcapon(tt)*G;
        Rx_M = Rx_M - loopGain*maxPeakValue_Pmusic(tt)*G;

    end

    Rxbf ( logical( eye( size(Rxbf,1)  ) ) ) = 0;
    Rx  ( logical( eye( size(Rx,1)   ) ) ) = 0;
    Rx_M( logical( eye( size(Rx_M,1) ) ) ) = 0;
    
    % Stop the iteration if the degraded CSM contains more information than
    % in the previous iteration
%     sumOfCSM = sum(sum(abs(Rx0)));
%     if sumOfCSM > sumOfDegradedCSM
%         break;
%     end
%     sumOfDegradedCSM = sumOfCSM;
    
end    
       
PbfC    = abs(PbfC);
PcaponC = abs(PcaponC + Pcapon_clean);
PmusicC = abs(PmusicC + Pmusic_clean);

end