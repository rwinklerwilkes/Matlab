#Forward/Backward Difference
function fd = forward(xZero,f,h)
	fd = (1/h)*(f(xZero)-f(xZero+h));
endfunction

#Three-Point Endpoint Formula
function en = end(xZero,f,h)
	en = (1/(2*h))*(-3*f(xZero)+4*f(xZero+h)-f(xZero+2*h));
	return;
endfunction

#Three-Point Midpoint Formula
function mi = mid(xZero,f,h)
	mi = (1/(2*h))*(f(xZero+h)-f(xZero-h));
	return;
endfunction

#Five-Point Endpoint Formula
function en = fend(xZero,f,h)
	en = (1/(12*h))*(-25*f(xZero)+48*f(xZero+h)-36*f(xZero+2*h)+16*f(xZero+3*h)-3*f(xZero+4*h));
	return;
endfunction

#Five-Point Midpoint Formula
function mi = fmid(xZero,f,h)
	mi = (1/(12*h))*(f(xZero-2*h)-8*f(xZero-h)+8*f(xZero+h)-f(xZero+2*h));
	return;
endfunction

#Composite Simpson's Rule
function s = simpson(f, a, b, n)
	h = (b-a)/n;
	aZero = f(a)+f(b);
	aOne = 0;
	aTwo = 0;
	i = 1;
	while(i<=n-1)
		X = a+i*h;
		if(mod(i,2)==0)
			aTwo = aTwo + f(X);
		else
			aOne = aOne + f(X);
		endif
		i++;
	endwhile
	s = (h/3)*(aZero + 2*aTwo + 4*aOne);
	return;
endfunction

#Composite Simpson's Rule (with tolerance)
function s = simpval(f, a, b, minTol,n)
	val = simpson(f,a,b,n);
	nextVal = simpson(f,a,b,n+2);
	n = n+2;
	while(abs(val-nextVal)>=minTol)
		n = n+2;
		val = nextVal;
		nextVal = simpson(f,a,b,n);
	endwhile
	disp("Nodes: "),disp(n);
	s = val;
	return;
endfunction

#Composite Trapezoidal Rule
function t = trapezoid(f,a,b,n)
	h = (b-a)/n;
	aZero = 0;
	aOne = 0;
	i = 0;
	while(i<=n-1)
		cur = a+i*h;
		next = a+(i+1)*h;
		aZero = aZero + h/2*(f(cur)+f(next));
		i++;
	endwhile
	t = aZero;
	return;
endfunction

#Composite Midpoint Rule
function m = midpoint(f,a,b,n)
	h = (b-a)/(n+2);
	i = 0;
	aZero = 0;
	while(i<(n/2+1))
		cur = a+h*(2*i+1);
		aZero = aZero+2*h*f(cur);
		i++;
	endwhile
	m = aZero;
	return;
endfunction

#Romberg Extrapolation
function r = romberg(f,a,b,n)
	rVals = zeros(n,n);
	i = 1;
	while(i<=n)
		rVals(i,1) = trapezoid(f,a,b,2**(i-1));
		i++;
	endwhile
	i = 2;
	while(i<=n)
		j = i;
		while(j<=n)
			rVals(j,i) = rVals(j,i-1)+(1/(4**(i-1)-1))*(rVals(j,i-1)-rVals(j-1,i-1));
			j++;
		endwhile
		i++;
	endwhile
	r = rVals;
	return;
endfunction

#Romberg Extrapolation
function r = rommin(f,a,b,maxIter,tolerance)
	n = 3;
	lastVal = romberg(f,a,b,n)(n,n);
	while(n<=maxIter)
		n++;
		value = romberg(f,a,b,n)(n,n);
		if(abs(lastVal-value)<tolerance)
			r = [value,n];
			return;
		endif
		lastVal = value;
	endwhile
	disp("Failed!");
	return;
endfunction

#Simpson's Method (for Use in Adaptive Quadrature)
function s = simp(f,a,b)
	mid = (b-a)/2;
	s = (mid/3)*(f(a)+4*f(a+mid)+f(b));
	return;
endfunction

#Adaptive Quadrature
function ad = adapt(f,a,b,tol)
	first = simp(f,a,(a+b)/2);
	second = simp(f,(a+b)/2,b);
	if(abs(simp(f,a,b)-first-second) < tol)
		ad = first+second;
	else
		ad = adapt(f,a,(a+b)/2,tol) + adapt(f,(a+b)/2,b,tol);
	endif
	return;
endfunction

#Simpson's Double Integral
function s = sdi(a,b,m,n,cx,dx,fxy)
	h = (b-a)/n;
	Jend = 0;
	Jeven = 0;
	Jodd = 0;
	i = 0;
	while(i<=n)
		x = a+i*h;
		HX = (dx(x)-cx(x))/m;
		Kend = fxy(x,cx(x)) + fxy(x,dx(x));
		Keven = 0;
		Kodd = 0;
		j = 1;
		while(j<=m-1)
			y = cx(x)+j*HX;
			Q = fxy(x,y);
			if(mod(j,2)==0)
				Keven = Keven + Q;
			else
				Kodd = Kodd + Q;
			endif
			j++;
		endwhile
		L = (HX/3)*(Kend + 2*Keven + 4*Kodd);
		if(i==0||i==n)
			Jend = Jend + L;
		elseif(mod(i,2)==0)
			Jeven = Jeven + L;
		else
			Jodd = Jodd + L;
		endif
		i++;
	endwhile
	s = (h/3)*(Jend + 2*Jeven + 4*Jodd);
	return;
endfunction

#Gaussian Double Integral
#Since I couldn't manage to figure out a way for Octave to generate these values, I hard-coded them in
#This restricts the largest m,n values to 5, however
function g = gdi(a,b,m,n,c,d,f)
	if(m>5||n>5)
		disp("Not enough roots in memory");
		return;
	endif
	roots = [0,0,0,0,0;.5773502692,-.5773502692,0,0,0;.7745966692,0,-.7745966692,0,0;.8611363116,.3399810436,-.3399810436,-.8611363116,0;.9061798459,.5384693101,0,-.5384693101,-.9061798459];
	coefficients = [0,0,0,0,0;1,1,0,0,0;.555555555,.888888888,.555555555,0,0;.3478548451,.6521451549,.6521451549,.3478548451,0;.2369268850,.4786286705,.56888888888,.4786286705,.2369268850];
	h1 = (b-a)/2;
	h2 = (b+a)/2;
	J = 0;
	i = 1;
	while(i<=m)
		JX = 0;
		x = h1*roots(m,i) + h2;
		d1 = d(x);
		c1 = c(x);
		k1 = (d1-c1)/2;
		k2 = (d1+c1)/2;
		j = 1;
		while(j<=n)
			y = k1*roots(n,j) + k2;
			Q = f(x,y);
			JX = JX + coefficients(n,j)*Q;
			j++;
		endwhile
		J = J + coefficients(m,i)*k1*JX;
		i++;
	endwhile
	g = h1*J;
	return;
endfunction

#Gaussian Triple Integral
#This has the same problem that the Gaussian double integral had
#Also, it doesn't work, I can't figure out why
function g = gti(a,b,m,n,p,c,d,alpha,beta,f)
	if(m>5||n>5||p>5)
		disp("Not enough roots in memory");
		return;
	endif
	roots = [0,0,0,0,0;.5773502692,-.5773502692,0,0,0;.7745966692,0,-.7745966692,0,0;.8611363116,.3399810436,-.3399810436,-.8611363116,0;.9061798459,.5384693101,0,-.5384693101,-.9061798459];
	coefficients = [0,0,0,0,0;1,1,0,0,0;.555555555,.888888888,.555555555,0,0;.3478548451,.6521451549,.6521451549,.3478548451,0;.2369268850,.4786286705,.56888888889,.4786286705,.2369268850];
	h1 = (b-a)/2;
	h2 = (b+a)/2;
	J = 0;
	i = 1;
	while(i<=m)
		JX = 0;
		x = h1*roots(m,i) + h2;
		d1 = d(x);
		c1 = c(x);
		k1 = (d1-c1)/2;
		k2 = (d1+c1)/2;
		j = 1;
		while(j<=n)
			JY = 0;
			y = k1*roots(n,j) + k2;
			beta1 = beta(x,y);
			alpha1 = alpha(x,y);
			l1 = (beta1-alpha1)/2;
			l2 = (beta1+alpha1)/2;
			k = 1;
			while(k<=p)
				z = l1*roots(p,k)+l2;
				Q = f(x,y,z);
				JY = JY + coefficients(p,k)*Q;
				k++;
			endwhile
			JX = JX + coefficients(n,j)*l1*JY;
			j++;
		endwhile
		J = J + coefficients(m,i)*k1*JX;
		i++;
	endwhile
	J = h1*J;
	g = J;
	return;
endfunction