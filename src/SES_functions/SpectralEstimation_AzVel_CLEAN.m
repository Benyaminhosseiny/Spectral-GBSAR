% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
function [Pbf0,Pcapon0,Pmusic0,Pbf,Pcapon,Pmusic]=SpectralEstimation_AzVel_CLEAN(rc,w_r,lambda,ant_loc_az,time_vec,theta,vel,num_tar,mode,c_iter)
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

Pbf    =zeros([length(theta),length(vel)]);
Pcapon =zeros([length(theta),length(vel)]);
Pmusic =zeros([length(theta),length(vel)]);

Pbf_clean    =zeros([length(theta),length(vel)]);
Pcapon_clean =zeros([length(theta),length(vel)]);
Pmusic_clean =zeros([length(theta),length(vel)]);


%Initialise trimmed cross spectral matrix (CSM) by setting the diagonal to zero
Rx0(logical(eye(size(Rx0,1)))) = 0;
Rx(logical(eye(size(Rx,1)))) = 0;

%Initialise break criterion
sumOfCSM = sum(sum(abs(Rx0)));
sumOfDegradedCSM = sumOfCSM;

Rx_M = Rx0;
% Rx_M = Rx;
Rxbf = Rx;
% Rxbf = Rx0;

[eigenVec,eigenVal] = eig(Rx_M);
Vn = eigenVec(:,1:end-1); % Noise Subspace

maxIterations = c_iter; loopGain = 0.9; % CLEAN parameters
for cleanMapIterations = 1:maxIterations
    Rx_inv = inv(Rx);
    if cleanMapIterations>1
        [eigenVec,eigenVal] = eig(Rx_M);
        Vn = eigenVec(:,1:end-1);
    end
    
    for thet_ii = 1:length(theta)
        for vel_ii = 1:length(vel)
            % Steering (Sensing) vector:
            SS = URA_Az_Vel(theta(thet_ii),vel(vel_ii));
            SS = SS(:);
            
            % BF:
            PP = SS'*Rxbf*SS;
            Pbf(thet_ii,vel_ii) = PP;
            %Capon:
            PP = SS'*Rx_inv*SS;
%             PP = SS'/Rx*SS;
            Pcapon(thet_ii,vel_ii) = 1/PP;
            % MUSIC:
            PP = SS'*(Vn*Vn')*SS;
            Pmusic(thet_ii,vel_ii) = 1/ PP;
        end
    end
    
    if cleanMapIterations == 1
        Pbf0    = abs(Pbf);
        Pcapon0 = abs(Pcapon);
        Pmusic0 = abs(Pmusic);
    end
    %% CLEAN:
    % -------------------------------------------------------
    % 2. Find peak value and its position in dirty map
    %     [maxPeakValue, maxPeakIndx] = max( abs( Pbf(:) ) );
    %     [maxPeakValueIndx_thet, maxPeakValueIndx_A, maxPeakValueIndx_f] = ind2sub( size(Pbf), maxPeakIndx );
    
    
    % Find peaks BF
    [peaks, locations] = findpeaks( abs( Pbf(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pbf, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pbf = locations(sortedIndices);
    
    % Find peaks Capon
    [peaks, locations] = findpeaks( abs( Pcapon(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pcapon, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pcapon = locations(sortedIndices);
    
    % Find peaks MUSIC
    [peaks, locations] = findpeaks( abs( Pmusic(:) ) );
    % Sort values and locations from biggest to smallest
    [maxPeakValue_Pmusic, sortedIndices] = sort(peaks, 'descend');
    sortedLocations_Pmusic = locations(sortedIndices);
    
    % Clean beam with specified width and max value of 1
    PmaxCleanBeam = zeros(length(theta), length(vel));
    g = 0;
    for tt = 1:num_tar % Three scatterers
        [maxPeakValueIndx_thet(tt), maxPeakValueIndx_vel(tt)] = ind2sub( size(Pbf), sortedLocations_Pbf(tt) );
        
        % -------------------------------------------------------
        % 3. Calculate the CSM induced by the peak source
        
        % Steering vector to location of peak source
        g = 0+URA_Az_Vel( theta(maxPeakValueIndx_thet(tt)), vel(maxPeakValueIndx_vel(tt)) );
        g = g(:);
        % Cross spectral matrix induced by peak source in that direction (eq. 11)
        G = g*g';
        G(logical(eye(size(G,1)))) = 0;
        
        % -------------------------------------------------------
        % 4. New updated map with clean beam from peak source location
        % Clean beam with specified width and max value of 1
        
        PmaxCleanBeam(maxPeakValueIndx_thet(tt), maxPeakValueIndx_vel(tt)) = 1;
        
        % Update clean map with clean beam from peak source location
        Pbf_clean    = Pbf_clean    + loopGain*maxPeakValue_Pbf(tt)   *PmaxCleanBeam;
        Pcapon_clean = Pcapon_clean + loopGain*maxPeakValue_Pcapon(tt)*PmaxCleanBeam;
        Pmusic_clean = Pmusic_clean + loopGain*maxPeakValue_Pmusic(tt)*PmaxCleanBeam;
        
        % -------------------------------------------------------
        % 5. Calculate degraded cross spectral matrix
        % Basically removing the PSF from that location of the plot
        Rxbf  = Rxbf  - loopGain*maxPeakValue_Pbf (tt)*G;
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

Pbf    = abs(Pbf);
Pcapon = abs(Pcapon + Pcapon_clean);
Pmusic = abs(Pmusic + Pmusic_clean);


end