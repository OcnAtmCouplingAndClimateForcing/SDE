function s=AR1_model(time,f,gamma)
% function s=AR1_model(time,f,gamma)
% TIME and F are the time index and forcing timseries for the AR1 model.
% GAMMA is the 1/(decorrelation timescale) following the AR1 equation:
% S(t+1) = (1-GAMMA*DT)*S(t) + F(t)*DT
%
% Output of the model is structured array S
% S.time  --> time index
% S.sig   --> integration of F following the AR1 model eq.
%
% E. Di Lorenzo (edl@gatech.edu)

% make sure the forcing has zero mean
f=f-mean(f);

% initialize variables
T=length(f);
sig=zeros(T,1);
sig(1)=0;

for t=1:T
  sig(t+1) = sig(t) *(1-gamma) + f(t);
end

% remove trends
sig=detrend(sig);

% interpolate output on the original time grid
s.sig=(sig(2:end)+sig(1:end-1))/2;
s.time=time;

% normalize the output by STD
s.sig=s.sig/std(s.sig);

