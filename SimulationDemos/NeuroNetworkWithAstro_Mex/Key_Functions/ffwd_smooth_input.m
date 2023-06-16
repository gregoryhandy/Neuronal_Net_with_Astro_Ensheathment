%%
% This function creates the feedforward input the goes into the non-spatial
% model
%%
function [Ix1e, Ix1i, Ix2e, Ix2i] = ffwd_smooth_input(num_ffwd_inputs,...
    taux,dt,Nt,mxe,mxi,sigma_s)

% Kx needs to be the sqrt of the desired auto-covariance
% in the sense of convolutions, ie Kx*Kx=Ax...the math works!
Kx=sqrt((1/taux)*sqrt(2/pi))*exp(-(-6*taux:dt:6*taux).^2./((taux)^2));
if num_ffwd_inputs == 1
    white_noise = [zeros(10000,1); randn(Nt-10000,1)./sqrt(dt)];
else
    white_noise = randn(Nt,1)./sqrt(dt);
end
sig1=dt*conv(white_noise,Kx,'same'); % white noise convolved with Kx

Ix1e=mxe+sqrt(sigma_s)*sig1;
Ix1i=mxi+sqrt(sigma_s)*sig1;

if num_ffwd_inputs == 1
    sig2 = sig1; % same input for both populations
else
    % Independent realization for population 2
    white_noise = randn(Nt,1)./sqrt(dt);
    sig2=dt*conv(white_noise,Kx,'same'); % white noise convolved with Kx
end
Ix2e=mxe+sqrt(sigma_s)*sig2;
Ix2i=mxi+sqrt(sigma_s)*sig2;

end