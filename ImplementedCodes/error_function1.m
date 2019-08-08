function [out]=error_function1(X,x,y)
% y2=X(1)*x.^2+X(2)*x+X(3);
y2=X(1).*exp(X(2)*x)+X(3)*x;

out=sum((y-y2).^2);

