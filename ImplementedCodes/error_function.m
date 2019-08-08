function [out]=error_function(X,x,y)
y2=X(1)*x.^2+X(2)*x+X(3); % function with 3 parameters X1, X2 and X3
out=sum((y-y2).^2);


