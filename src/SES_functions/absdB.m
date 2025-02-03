function [output] = absdB(input)
% input(input==0)=1e-20;
output=10*log10(abs(input));
end

