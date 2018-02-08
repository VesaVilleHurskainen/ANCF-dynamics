function [ status ] = OutputFcnWrapper( t, y, flag, data )
%OUTPUTFCNWRAPPER Wrapper for OutputFcn for use with radau5.
%   Detailed explanation goes here
if strcmp(flag,'interp')
    % radau5 code only works as intended if output only happens for flag
    % 'itnerp'
    status = OutputFcn( t, y, flag, data );
else
    status=0;
end
end

