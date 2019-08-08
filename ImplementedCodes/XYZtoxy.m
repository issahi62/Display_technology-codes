% %  function to calculate x y z chromaticity coordinates 


function [output1]=XYZtoxy(XYZ)

% calculate x and y and z
x=XYZ(:,1)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3)); 
y=XYZ(:,2)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3));
z=XYZ(:,3)./(XYZ(:,1)+XYZ(:,2)+XYZ(:,3));

output1(:,1)=x;
output1(:,2)=y;
output1(:,3)=z;