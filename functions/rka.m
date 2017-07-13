function [xSmall, t, dt] = rka(x,t,dt,err,derivsRK,param)
% Adaptive Runge-Kutta routine
% Inputs
%   x          Current value of the dependent variable
%   t          Independent variable (usually time)
%   dt        Step size (usually time step)
%   err        Desired fractional local truncation error
%   derivsRK   Right hand side of the ODE; derivsRK is the
%              name of the function which returns dx/dt
%              Calling format derivsRK(x,t,param).
%   param      Extra parameters passed to derivsRK
% Outputs
%   xSmall     New value of the dependent variable
%   t          New value of the independent variable
%   dt        Suggested step size for next call to rka

%* Set initial variables
tSave = t;  xSave = x;    % Save initial values
safe1 = .9;  safe2 = 4.;  % Safety factors

%* Loop over maximum number of attempts to satisfy error bound
maxTry = 100;  
for iTry=1:maxTry
	
  %* Take the two small time steps
  half_dt = 0.5 * dt;
  xTemp = rk4(xSave,tSave,half_dt,derivsRK,param);
  t = tSave + half_dt;
  xSmall = rk4(xTemp,t,half_dt,derivsRK,param);
  
  %* Take the single big time step
  t = tSave + dt;
  xBig = rk4(xSave,tSave,dt,derivsRK,param);
  
  %* Compute the estimated truncation error
  scale = err * (abs(xSmall) + abs(xBig))/2.;
  xDiff = xSmall - xBig;
  errorRatio = max( abs(xDiff)./(scale + eps) );
  
  %* Estimate new dt value (including safety factors)
  dt_old = dt;
  dt = safe1*dt_old*errorRatio^(-0.20);
  dt = max(dt,dt_old/safe2);
  dt = min(dt,safe2*dt_old);
  
  %* If error is acceptable, return computed values
  if (errorRatio < 1)  return;  end 
end

%* Issue error message if error bound never satisfied
error('ERROR: Adaptive Runge-Kutta routine failed');
return;
  
