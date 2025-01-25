% Citation:
% Hosseiny, Benyamin, Jalal Amini, and Hossein Aghababaei. 
% "Spectral estimation model for linear displacement and vibration monitoring with GBSAR system." 
% Mechanical Systems and Signal Processing 208 (2024): 110916.
% https://doi.org/10.1016/j.ymssp.2023.110916
% Numerical simulation of SAR/MIMO signal
%   * Latest Update: 2023-09-12 *

%% SIGNAL
t = linspace(0,T,Nr)';
beam_az = rad2deg(lambda/m_az/d_az);
beam_el = rad2deg(lambda/m_el/d_el);
signal_MIMO_TS = zeros(Nr,Na);


for tar_ii = 1:num_tar
    %% Polar:
    if strcmp('POLAR',upper(signal_model))
        R_tar_ii = R_tar(tar_ii)+dR_tar_total_disp(:,tar_ii)';
        tau = 2*R_tar_ii/c;
        signal_MIMO_TS = signal_MIMO_TS+ ...
            exp( 1i*(...
            2*pi*(fc*tau + cr*t*tau - cr*(tau.^2)/2) ...
            + m_az*2*pi*X_rad*sind(theta_tar(tar_ii))/lambda ...
            + (m_az-1)*4*pi*cr*t.*X_rad*sind(theta_tar(tar_ii))/c ...       % Azimuth
            + m_el*2*pi*(Z_rad*sind(ph_tar(tar_ii))/lambda)...              % Elevation
            ) );
    end
    %% Cartesian:
    if strcmp('CARTESIAN',upper(signal_model))
        R_tar_ii = sqrt( (X_tar(tar_ii)-X_rad).^2+(Y_tar(tar_ii)-Y_rad).^2+(Z_tar(tar_ii)-Z_rad).^2 );
        tau = (m_az*R_tar_ii-2*dR_tar_total_disp(:,tar_ii)')/c;
        signal_MIMO_TS = signal_MIMO_TS + exp( -1i*2*pi*( fc*tau + cr*t*tau - cr*(tau.^2)/2 ) );
               
        
    end
end
signal_MIMO_TS = awgn( signal_MIMO_TS,snr );