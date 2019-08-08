function XYZ=DigittoXYZ_GOG(r,g,b,offset,test,CMFs)

r=r-ones(18,1)*offset;
g=g-ones(18,1)*offset;
b=b-ones(18,1)*offset;


R=sum(r,2)./sum(r(18,:));
G=sum(g,2)./sum(g(18,:));
B=sum(b,2)./sum(b(18,:));

% Check if initial values are well selected

% because we have 8 bits 0-255
ROy=fitting_GOG(0:15:255,R',[1,.1,1]);
GOy=fitting_GOG(0:15:255,G',[1,.1,1]);
BOy=fitting_GOG(0:15:255,B',[1,.55,1]);

% Check if initial values are well selected


% xx=0:15:255
% tmp=(ROy(1)*(xx/(2^8-1))+ROy(2)).^ROy(3);
% 
% figure()
% plot(0:15:255,R,'o',xx,tmp,'r')
% title('Red channel')
% 
% tmp=(GOy(1)*(xx/(2^8-1))+GOy(2)).^GOy(3);
% 
% figure()
% plot(0:15:255,G,'o',xx,tmp,'g')
% title('Green channel')
% 
% tmp=(BOy(1)*(xx/(2^8-1))+BOy(2)).^BOy(3);
% 
% figure()
% plot(0:15:255,B,'o',xx,tmp,'b')
% title('Blue channel')


%Estimate RGB value from the digit value dr dg db using GOG model
R_dr=(ROy(1)*(test(1)/(2^8-1))+ROy(2)).^ROy(3);
G_dg=(GOy(1)*(test(2)/(2^8-1))+GOy(2)).^GOy(3);
B_db=(BOy(1)*(test(3)/(2^8-1))+BOy(2)).^BOy(3);

% Estimate XYZ values for maximum R G B
R_max=CalculateXYZ(r(18,:),CMFs);
G_max=CalculateXYZ(g(18,:),CMFs);
B_max=CalculateXYZ(b(18,:),CMFs);

% Estimate XYZ values for offset
XYZoffset=CalculateXYZ(offset,CMFs);

% Build Matrix
Matrix_Max=[R_max;G_max;B_max]';
Matrix_Digit=[R_dr,G_dg,B_db]';
Matrix_Offset=XYZoffset';

% Estimate XYZ
XYZ=Matrix_Max*Matrix_Digit+Matrix_Offset;

