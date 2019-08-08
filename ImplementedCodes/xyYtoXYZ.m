function [output]= xyYtoXYZ(xyY)
X=xyY(:,3)./xyY(:,2).*xyY(:,1);
z=ones(size(xyY(:,1)))-xyY(:,1)-xyY(:,2);
Z=xyY(:,3)./xyY(:,2).*z;

output=[X xyY(:,3) Z];