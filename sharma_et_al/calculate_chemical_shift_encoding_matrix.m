function A = calculate_chemical_shift_encoding_matrix(FieldStrength, ppm, TE)

% Author:  Samir Sharma
% Created: November 2011

GAMMA = 42.58; %in MHz/T
for i = 1:length(ppm)
  df = GAMMA * FieldStrength * ppm(i);
  A(:,i) = exp(1i*2*pi*TE(:)*df); % DH*: removed "-" sign from exponential to conform with toolbox
%  A(:,aa) = exp(-1i*2*pi*data.TE(:)*df) * algoParams.species(aa).relAmps(:);
end