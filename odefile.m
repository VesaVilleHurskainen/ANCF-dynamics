function [ ydot ] = odefile( t,y,options,data )
%ODEFILE wrapper for function eom_Dynam for use with radau5 code
[ydot, ~, ~] = eom_Dynam(t,y,data);

end

