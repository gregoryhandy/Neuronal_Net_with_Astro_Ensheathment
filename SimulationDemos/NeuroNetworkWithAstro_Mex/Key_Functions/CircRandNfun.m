%%
% This function generates n random numbers that are integers 
% between rand_min and rand_max.
% The distribution of these random variables is like a Gaussian 
% distribution with mean mu and std sigma, but rounded to the nearest 
% integer and wrapped around the interval [min,max]. 
%%%
function [ z ] = CircRandNfun(mu, sigma, rand_min, rand_max, n ) 

u1 = round(mu + sigma*randn(n,1))-rand_min;
z = mod(u1,rand_max-rand_min+1)+rand_min;

end

% Non-vector version of Rosenbaum's code
% for i=1:2:(n-1)
%     u1=rand();
%     u2=rand();
%      
%     lower_num_1 = round(sigma*sqrt(-2*log(u1))*cos(2*pi*u2)+mu)-rand_min;
%     lower_num_2 = round(sigma*sqrt(-2*log(u1))*sin(2*pi*u2)+mu)-rand_min;
%     z(i) = mod(lower_num_1,rand_max-rand_min+1)+rand_min;
%     z(i+1) = mod(lower_num_2,rand_max-rand_min+1)+rand_min;   
% end
% 
% 
% if(i==(n-2))
%     u1=rand();
%     u2=rand();
%     lower_num_1 = round(sigma*sqrt(-2*log(u1))*cos(2*pi*u2)+mu)-rand_min;
%     z(i+2)=mod(lower_num_1,rand_max-rand_min+1)+rand_min;
% end