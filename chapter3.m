#Chapter 3 - Numerical Analysis
#Rich Winkler

#Neville's Method
function n = neville(xVals, qVals, xPoint)
	iter = 1;
	n = zeros(columns(xVals),columns(xVals));
	while(iter<=columns(xVals))
		n(iter,1) = qVals(iter);
		iter++;
	endwhile
	iter = 2;
	while(iter<=columns(xVals))
		jter = 2;
		while(jter<=iter)
			newVal = (1/(xVals(iter)-xVals(iter-jter+1)))*((xPoint-xVals(iter-jter+1))*n(iter,jter-1)-(xPoint-xVals(iter))*n(iter-1,jter-1));
			n(iter,jter)=newVal;
			jter++;
		endwhile
		iter++;
	endwhile
	return;
endfunction

#Divided-Difference Formula
function d = fdiff(xVals, fVals)
	d = zeros(columns(xVals),columns(xVals));
	i = 1;
	while(i<=columns(xVals))
		d(i,1)=fVals(i);
		i++;
	endwhile
	i = 2;
	while(i<=columns(xVals))
		j = 2;
		while(j <= i)
			d(i,j) = (1/(xVals(i) - (xVals(i-j+1))))*(d(i,j-1)-d(i-1,j-1));
			j++;
		endwhile
		i++;
	endwhile
	return;
endfunction

#Backward Divided-Difference Formula
function b = bdiff(xVals,fVals)
	b = zeros(columns(xVals),columns(xVals));
	i = 1;
	while(i<=columns(xVals))
		b(i,1)=fVals(columns(xVals)-i+1);
		i++;
	endwhile
	i = 2;
	while(i<=columns(xVals))
		j = 2;
		while(j <= i)
			b(i,j) = (1/(xVals(i-j+1) - (xVals(i))))*(b(i,j-1)-b(i-1,j-1));
			j++;
		endwhile
		i++;
	endwhile
	return;
endfunction

#Hermite Polynomial
function h = hermite(xVals, fVals, primeVals)
	iter = 1;
	z = linspace(0,0,2*columns(xVals));
	Q = zeros(2*columns(xVals));
	while(iter<=columns(xVals))
		z(2*iter-1) = xVals(iter);
		z(2*iter) = xVals(iter);
		Q(2*iter-1,1) = fVals(iter);
		Q(2*iter,1) = fVals(iter);
		Q(2*iter,2) = primeVals(iter);
		if(iter>1)
			Q(2*iter-1,2) = (Q(2*iter-1,1)-Q(2*iter-2,1))/(z(2*iter-1)-z(2*iter-2));
		endif
		iter++;
	endwhile
	i = 3;
	while(i<=2*columns(xVals))
		j = 3;
		while(j<=i)
			Q(i,j) = (Q(i,j-1)-Q(i-1,j-1))/(z(i)-z(i-j+1));
			j++;
		endwhile
		i++;
	endwhile
	h = Q;
	return;
endfunction

#Cubic Spline Interpolation - Natural
function s = cspline(xVals,fVals)
	iter = 1;
	n = columns(xVals);
	h = linspace(0,0,n);
	while(iter<=(n-1))
		h(iter) = xVals(iter+1)-xVals(iter);
		iter++;
	endwhile
	iter = 2;
	alpha = linspace(0,0,n);
	while(iter<=n-1)
		alpha(iter) = (3/(h(iter)))*(fVals(iter+1)-fVals(iter))-(3/h(iter-1))*(fVals(iter)-fVals(iter-1));
		iter++;
	endwhile
	l = linspace(0,0,n);
	mu = linspace(0,0,n);
	z = linspace(0,0,n);
	iter = 2;
	l(1) = 1;
	while(iter<=n-1)
		l(iter) = 2*(xVals(iter+1)-xVals(iter-1))-h(iter-1)*mu(iter-1);
		mu(iter) = (h(iter)/l(iter));
		z(iter) = (alpha(iter)-h(iter-1)*z(iter-1))/l(iter);
		iter++;
	endwhile
	l(n) = 1;
	z(n) = 0;
	d = linspace(0,0,n);
	b = linspace(0,0,n);
	c = linspace(0,0,n);
	c(n) = 0;
	iter = n-1;
	while(iter>=1)
		c(iter) = z(iter)-mu(iter)*c(iter+1);
		b(iter) = (fVals(iter+1)-fVals(iter))/h(iter)-(h(iter)*(c(iter+1)+2*c(iter)))/3;
		d(iter) = (c(iter+1)-c(iter))/(3*h(iter));
		iter--;
	endwhile
	A = zeros(n,4);
	A(:,1) = fVals;
	A(:,2) = b;
	A(:,3) = c;
	A(:,4) = d;
	s = A;
	return;
endfunction

#CCubic Spline Interpolation - Clamped
function s = clspline(xVals,fVals,fPrimeZero,fPrimeN)
	iter = 1;
	n = columns(xVals);
	h = linspace(0,0,n);
	while(iter<=(n-1))
		h(iter) = xVals(iter+1)-xVals(iter);
		iter++;
	endwhile
	iter = 2;
	alpha = linspace(0,0,n);
	alpha(1) = 3*(fVals(2) - fVals(1))/h(1)-3*fPrimeZero;
	alpha(n) = 3*fPrimeN-3*(fVals(n)-fVals(n-1))/h(n-1);
	while(iter<=n-1)
		alpha(iter) = (3/(h(iter)))*(fVals(iter+1)-fVals(iter))-(3/h(iter-1))*(fVals(iter)-fVals(iter-1));
		iter++;
	endwhile
	l = linspace(0,0,n);
	mu = linspace(0,0,n);
	z = linspace(0,0,n);
	iter = 2;
	l(1) = 2*h(1);
	mu(1) = .5;
	z(1) = alpha(1)/l(1);
	while(iter<=n-1)
		l(iter) = 2*(xVals(iter+1)-xVals(iter-1))-h(iter-1)*mu(iter-1);
		mu(iter) = (h(iter)/l(iter));
		z(iter) = (alpha(iter)-h(iter-1)*z(iter-1))/l(iter);
		iter++;
	endwhile
	l(n) = h(n-1)*(2-mu(n-1));
	z(n) = (alpha(n)-h(n-1)*z(n-1))/l(n);
	d = linspace(0,0,n);
	b = linspace(0,0,n);
	c = linspace(0,0,n);
	c(n) = z(n);
	iter = n-1;
	while(iter>=1)
		c(iter) = z(iter)-mu(iter)*c(iter+1);
		b(iter) = (fVals(iter+1)-fVals(iter))/h(iter)-(h(iter)*(c(iter+1)+2*c(iter)))/3;
		d(iter) = (c(iter+1)-c(iter))/(3*h(iter));
		iter--;
	endwhile
	A = zeros(n,4);
	A(:,1) = fVals;
	A(:,2) = b;
	A(:,3) = c;
	A(:,4) = d;
	s = A;
	return;
endfunction

function b = bezier(endpoint, leftGuides, rightGuides)
	n = rows(endpoint);
	a = zeros(n-1,4);
	c = zeros(n-1,4);
	iter = 1;
	while(iter<=n-1)
		a(iter,1) = endpoint(iter,1);
		c(iter,1) = endpoint(iter,2);
		a(iter,2) = 3*(leftGuides(iter,1)-endpoint(iter,1));
		c(iter,2) = 3*(leftGuides(iter,2)-endpoint(iter,2));
		a(iter,3) = 3*(endpoint(iter,1)+rightGuides(iter,1)-2*leftGuides(iter,1));
		c(iter,3) = 3*(endpoint(iter,2)+rightGuides(iter,2)-2*leftGuides(iter,1));
		a(iter,4) = endpoint(iter+1,1)-endpoint(iter,1)+3*leftGuides(iter,1)-3*rightGuides(iter,1);
		c(iter,4) = endpoint(iter+1,2)-endpoint(iter,2)+3*leftGuides(iter,2)-3*rightGuides(iter,2);
		iter++;
	endwhile
	b = zeros(n-1,8);
	b(:,1:4) = a;
	b(:,5:8) = c;
	return;
endfunction