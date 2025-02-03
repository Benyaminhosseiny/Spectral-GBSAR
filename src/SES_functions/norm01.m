function [output] = norm01(input)
output=(input-nanmin(input,[],'all'))./(nanmax(input,[],'all')-nanmin(input,[],'all'));
output=output+1e-10;
end