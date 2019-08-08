%  Calculate Tristimulus value
function [output] = CalculateXYZ(spectra,xyz_curves)
k = 683;
output = k*(spectra*xyz_curves);  % [X Y Z] as output
end
