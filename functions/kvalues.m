function k = kvalues(N)

% Returns the values of k used in FFT and IFFT.
% Why they look like this is a mystery and MATLAB is shit

if mod(N,2) == 0
	k = [0:N/2-1 0 -N/2+1:-1];
else
	error('N /HAS/ to be even, otherwise it wont work, because MATLAB is shit')
end
	