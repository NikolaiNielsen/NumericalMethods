function dudt = spectral(uin,t,p)

	k = p{1};
	D = p{2};
	tau = p{3};
	C0 = p{4};
	
	dudt = -(D*k.^2+1/tau).*uin+C0;

end