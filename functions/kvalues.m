function k = kvalues(N)

% Returns the values of k used in FFT and IFFT.
% Why they look like this is a mystery and MATLAB is shit

k = [0:N/2-1 0 -N/2+1:-1];