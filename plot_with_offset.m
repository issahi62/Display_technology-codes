function plot_with_offset(dat,xyz31_1nm)
hold on
c=xyz31_1nm(21:2:421,2:4);
% plot each patch from the ramp, with a offset !!!
XYZ=CalculateXYZ(dat{4},c);
xy=XYZtoxy(XYZ);
plot(xy(:,1),xy(:,2),'mo'); 
legend('locus','sRGB','aRGB','screen', 'offset');
% You can see that chromaticity corresponding 
%to a small digital value is quite far away from the primary.
%It shows how big influence has offset signal to the overall chromaticity
hold off
end

