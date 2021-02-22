function noise = rician_noise(len, sd)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
s = sd;
v = 0;
noise = ricernd(v*ones(1,len), s);
end

