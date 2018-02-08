function [a,ad,add] = motion_case10(t)

A1 = 1;
T1 = 5;
omega = (1700/60)/(2188e-3/2);
% omega = 10000/60*2*pi;

if t<=T1
    a = A1*omega*(pi*t - T1*sin((pi*t)/T1))/(2*pi);
	ad = A1*omega*(1-cos(pi*t/T1))/2;
	add = (A1*omega*pi*sin((pi*t)/T1))/(2*T1);
else
    a = A1*omega*t-A1*omega*T1/2;
	ad = A1*omega;
	add = 0;
end 