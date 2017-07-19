function a = randmat(n,p)

% Creates an nxn matrix where approx p of the elements are 1
a = zeros(n);
a(rand(n)>=p) = 1;
a = xor(a,a');