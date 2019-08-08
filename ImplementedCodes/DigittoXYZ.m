function XYZ=DigittoXYZ(r,g,b,offset,test,CMFs)

% plot(380:2:780, r)
% CMFs=c

r=r-ones(18,1)*offset;
g=g-ones(18,1)*offset;
b=b-ones(18,1)*offset;


R=sum(r,2)./sum(r(18,:));
G=sum(g,2)./sum(g(18,:));
B=sum(b,2)./sum(b(18,:));

% plot(xx,R)

% dr=225;
xx=0:15:255;
R_dr=interp1(xx,R,test(1));
G_dg=interp1(xx,G,test(2));
B_db=interp1(xx,B,test(3));


R_max=CalculateXYZ(r(18,:),CMFs);
G_max=CalculateXYZ(g(18,:),CMFs);
B_max=CalculateXYZ(b(18,:),CMFs);

XYZoffset=CalculateXYZ(offset,CMFs);

Matrix_Max=[R_max;G_max;B_max]';
Matrix_Digit=[R_dr,G_dg,B_db]';
Matrix_Offset=XYZoffset';


XYZ=Matrix_Max*Matrix_Digit+Matrix_Offset;

