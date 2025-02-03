function dR_tri = displacement_model_triangle_TS(amp,freq,prt,time,fi_shift)
% Generate triangular displacement function by summation of HARMONICS.
% more: https://en.wikipedia.org/wiki/Triangle_wave
% amp : displacement amplitude (m)
% freq: triangle's fundamental frequency
% prt : data acquisition (sampling) interval
% time: data acquisition (sampling) total time
% fi_shift: starting phase shift

if nargin<5
    fi_shift = zeros(size(amp));
end
amp = amp(:)';
freq = freq(:)';
fi_shift = fi_shift(:)';
%

Nh = 3;2;1;7;         % number of harmonics [the higher the Nh the better approximation of tiangular signal]
ii = (0:Nh-1)'; % harmonic label
nn = 2*ii+1;    % harmonic mode number
dR_tri=[];
time_vec=prt*(0:round(time/prt));
for A_f_s_ii = [amp;freq;fi_shift]
    amp_i  = A_f_s_ii(1);
    freq_i = A_f_s_ii(2);
    fi_i   = A_f_s_ii(3);
    dR_tri = [ dR_tri, ...
                amp_i * sum( 8/(pi^2)*((-1).^ii).*(nn.^-2).*sin(2*pi*freq_i*nn*time_vec+fi_i+0*pi), 1 )' ];
end



% p: period (s)
% samples = round(t/prt);
% for ii = 0:samples
%     time_ii=prt*ii;
% %     dR(ii+1,:) = 2*amp.*( time_ii/p-floor(time_ii/p+0.5) );
%     dR_tri(ii+1,:) = amp.*( -2./p*abs(mod(time_ii,p)-p/2)+1 );
% end

end

