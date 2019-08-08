
function [ABC]=fitting(x,y,guess)

ABC=fminsearch(@error_function,guess,optimset('display','on'),x,y);  % find the minimizer for function error_function, starting points are fiven by guess variable