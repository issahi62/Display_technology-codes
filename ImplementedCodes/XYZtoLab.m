%% 12. CIE 1976 Lab space and color difference formula
% %colorconv.m 

function Lab=XYZtoLab(XYZ,XYZn)
X0=XYZn(1); Y0=XYZn(2); Z0=XYZn(3);

for n=1:size(XYZ,1)
FX=XYZ(n,1)/X0;
if FX > 0.008856
XX(1,n) = FX.^(1/3);
else
XX(1,n)=7.787*FX + 16/116;
end
FY=XYZ(n,2)/Y0;
if FY > 0.008856
YY(1,n) = FY.^(1/3);
else
YY(1,n)=7.787*FY + 16/116;
end
FZ=XYZ(n,3)/Z0;
if FZ > 0.008856
ZZ(1,n) = FZ.^(1/3);
else
ZZ(1,n)=7.787*FZ + 16/116;
end
end
Lt=116*YY-16;
at=500*(XX-YY);
bt=200*(YY-ZZ);
Lab=[Lt at bt];