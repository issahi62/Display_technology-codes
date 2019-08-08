% *************************
%Exercise 10 -color gamut codes for all the display in one graph.
%**************************
load( 'IN_LABORATORY_2019.mat', 'CRT', 'Dell24', 'Dell_Konica', 'EIZO_2012', 'Projector_DLP_BENQ', 'xyz31_1nm');

%******************
%Codes for plottion color gamut
%******************
data = eye(201);
data(202,:) = data(1,:);
cmf = xyz31_1nm(21:2:421,2:4);
XYZ1 = CalculateXYZ(data, cmf);
xy1 = XYZtoxy(XYZ1);
plot(xy1(:,1), xy1(:,2), 'k')
axis equal
grid on
hold on
XYZ2=CalculateXYZ(CRT{4}([18 36 54 18] ,:),cmf);
xy2 = XYZtoxy(XYZ2);
plot(xy2(:,1), xy2(:,2), 'b')
XYZ3=CalculateXYZ(Dell24{4}([18 36 54 18] ,:),cmf);
xy3 = XYZtoxy(XYZ3);
plot(xy3(:,1), xy3(:,2), 'r')
XYZ4=CalculateXYZ(Dell_Konica{4}([18 36 54 18] ,:),cmf);
xy4 = XYZtoxy(XYZ4);
plot(xy4(:,1), xy4(:,2), 'y')
XYZ5=CalculateXYZ(EIZO_2012{4}([18 36 54 18] ,:),cmf);
xy5 = XYZtoxy(XYZ5);
plot(xy5(:,1), xy5(:,2), 'g')
XYZ6=CalculateXYZ(Projector_DLP_BENQ{4}([18 36 54 18] ,:),cmf);
xy6 = XYZtoxy(XYZ6);
plot(xy6(:,1), xy6(:,2))
legend ('Data','CRT','Dell24','Dell Konica', 'EIZ0', 'Projector'); 
xlabel('x');
ylabel('y');
title('Chromaticity Diagram');