
function [ABC]=fitting1(x,y,guess)

ABC=fminsearch(@error_function1,guess,optimset('display','off'),x,y);