function plot_without_offset(dat,xyz31_1nm)
hold on
c=xyz31_1nm(21:2:421,2:4);
offsetRepmat=repmat(dat{4}(1,:),4,1); % Create matrix with  4 times Repeated offset for future calculations
XYZ = CalculateXYZ(dat{4}([18 36 54 18],:)-offsetRepmat,c); % Calculate Device's Gamut  (Now Offset is subtracted !!) 
xyY2 = XYZtoxy(XYZ);
plot(xyY2(:,1),xyY2(:,2),'g-o');

WithoutOffset = zeros(54,201); % create matrix full of '0'

for i=1:54
    WithoutOffset(i,:) = dat{4}(i,:) - dat{4}(1,:);  % Subtract offset from the each measurement
end
WithoutOffset(WithoutOffset<0)=0; % negative values assign to 0

% Chromacity constancy for LED MONITOR
XYZ = CalculateXYZ(WithoutOffset,xyz31_1nm(21:2:421,2:4));
xyY3 = XYZtoxy(XYZ);
plot(xyY3(:,1),xyY3(:,2),'bo');
legend('locus','sRGB','aRGB','screen', 'withoutoffset');
hold off
end