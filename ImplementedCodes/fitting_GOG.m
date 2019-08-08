
function [ABC]=fitting_GOG(x,y,guess)

ABC=fminsearch(@error_functionGOG,guess,optimset('display','off'),x,y);  % Function explaine in the Exercise 3