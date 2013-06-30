#Chapter 2 Source Code
#Rich Winkler

#Bisection method
function bisect = approx2( p , begl , endl , minTol ,maxIter )
	iter = 1;
	while(iter<=maxIter)
		mid = begl+(endl-begl)/2;
		if(p(mid) == 0 || ((endl-begl)/2<minTol))
			bisect = mid;
			return;
		endif
		iter++;
		if(p(begl)*p(mid)>0)
			begl = mid;
		else
			endl = mid;
		endif
	endwhile
endfunction

#Fixed-point iteration
function fix = fpoint( f, pzero, minTol, maxIter)
	iter = 1;
	while(iter<=maxIter)
		fix = f(pzero);
		if(abs(fix-pzero)<minTol)
			disp(iter);
			return;
		endif
		iter++;
		pzero = fix;
	endwhile
	disp("Failed");
	return;
endfunction

#Newton's method
function n = newton( f, fPrime, pZero, minTol, maxIter)
	iter = 1;
	while(iter<=maxIter)
		n = pZero - (f(pZero)/fPrime(pZero));
		if(abs(n-pZero)<minTol)
			disp("Iterations: "),disp(iter);
			return;
		endif
		pZero = n;
		iter++;
	endwhile
	disp("Failed");
endfunction

#Secant method
function s = secmet( f, pZero, pOne, minTol, maxIter)
	iter = 1;
	qZero = f(pZero);
	qOne = f(pOne);
	while(iter<=maxIter)
		s= pOne - ((qOne*(pOne-pZero))/(qOne-qZero));
		if(abs(s-pOne)<minTol)
			disp("Iterations: "), disp(iter);
			return;
		endif
		pZero = pOne;
		qZero = qOne;
		pOne = s;
		qOne = f(s);
		iter++;
	endwhile
	disp("Failed");
endfunction

#Method of False Position
function fal = falsepos(f, pZero, pOne, minTol, maxIter)
	iter = 1;
	qZero = f(pZero);
	qOne = f(pOne);
	while(iter<=maxIter)
		fal= pOne - ((qOne*(pOne-pZero))/(qOne-qZero));
		if(abs(fal-pOne)<minTol)
			disp("Iterations: "), disp(iter);
			return;
		endif
		q = f(fal);
		if(q*qOne<0)
			pZero = pOne;
			qZero = qOne;
		endif
		pOne = fal;
		qOne = q;
		iter++;
	endwhile
	disp("Failed");
endfunction

#Aitken's Method
#Roughly twice as slow as Steffensen's on x^3+4x^2-10=0
function a = aitken(f, pZero, minTol, maxIter)
	iter = 1;
	while(iter<=maxIter)
		pOne = f(pZero);
		pTwo = f(pOne);
		a = pZero-(((pOne-pZero)**2)/(pTwo-2*pOne+pZero));
		if(abs(a-pZero)<minTol)
			disp("Iterations: "),disp(iter)
			return;
		endif
		pZero = pOne;
		iter++;
	endwhile
endfunction

#Steffensen's Method
function s = steffensen(f, pZero, minTol, maxIter)
	iter = 1;
	while(iter<=maxIter)
		pOne = f(pZero);
		pTwo = f(pOne);
		s = pZero-(((pOne-pZero)**2)/(pTwo-2*pOne+pZero));
		if(abs(s-pZero)<minTol)
			disp("Iterations: "),disp(iter)
			return;
		endif
		pZero = s;
		iter++;
	endwhile
	disp("Failed");
endfunction

#Horner's Method
#This is a modified version, which allows you to increase the number of iterations for more accuracy
function h = horner(poly, initVal, numIter)
	iter = 1;
	while(iter<=numIter)
		y = poly(1);
		z = poly(1);
		for i = 2:(columns(poly)-1)
			y = initVal*y+poly(i);
			z = initVal*z+y;
		endfor
		y = initVal*y+poly(columns(poly));
		h = initVal-(y/z);
		initVal = h;
		iter++;
	endwhile
	return;
endfunction

#Muller's Method
function m = muller(f, pZero,pOne,pTwo, minTol, maxIter)
	hOne = pOne-pZero;
	hTwo = pTwo-pOne;
	dOne = (f(pOne)-f(pZero))/hOne;
	dTwo = (f(pTwo)-f(pOne))/hTwo;
	d = (dTwo-dOne)/(hTwo+hOne);
	iter = 3;
	while(iter<=maxIter)
		b = dTwo+hTwo*d;
		D = (b**2-4*f(pTwo)*d)**.5;
		if(abs(b-D)<abs(b+D))
			E = b+D;
		else
			E = b-D;
		endif
		h = (-2*f(pTwo))/E;
		p = pTwo+h;
		if(abs(h)<minTol)
			disp("Iterations: "), disp(iter);
			m = p;
			return;
		endif
		pZero = pOne;
		pOne = pTwo;
		pTwo = p;
		hOne = pOne-pZero;
		hTwo = pTwo-pOne;
		dOne = (f(pOne)-f(pZero))/hOne;
		dTwo = (f(pTwo)-f(pOne))/hTwo;
		d = (dTwo-dOne)/(hTwo+hOne);
		iter++;
	endwhile
	disp("Failed");
endfunction