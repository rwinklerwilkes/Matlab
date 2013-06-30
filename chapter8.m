function l = lsa(x,y)
	n = columns(x);
	if(n!=columns(y))
		disp("x must have same dimension as f(x)");
		return;
	endif
	sumx = 0;
	sumy = 0;
	sumxy = 0;
	sumxsq =0;
	for i = 1:n
		sumx = sumx + x(i);
		sumy = sumy + y(i);
		sumxy = sumxy + x(i)*y(i);
		sumxsq = sumxsq + x(i)**2;
	endfor
	a = zeros(1,2);
	a(1) = (sumxsq*sumy-sumxy*sumx)/(n*sumxsq-sumx**2);
	a(2) = (n*sumxy-sumx*sumy)/(n*sumxsq-sumx**2);
	l = a;
	return;
endfunction

