function [out]=error_functionGOG(X,x,y)
% y2=X(1)*x.^2+X(2)*x+X(3); ex.1
% y2=X(1).*exp(X(2)*x)+X(3)*x; ex.2

N=8;
y2=(X(1)*(x/(2^N-1))+X(2)).^X(3); % ex.4   


out=sum((y-y2).^2);

