%Power Method
function [p,eigval] = pm(A,x,TOL,maxIter)
	k = 1;
	p = 1;
	n = columns(A);
	ixp = 0;
	ix = infnorm(x);
	while(p<=n)
		if(abs(x(p))==ix)
			ixp = p;
			break;
		endif
		p++;
	endwhile
	x = x/ix;
	while(k<=maxIter)
		y = A*x;
		mu = y(ixp);
		p = 1;
		iy = infnorm(y);
		iyp = 0;
		while(p<=n)
			if(abs(y(p))==iy)
				iyp = p;
				break;
			endif
			p++;
		endwhile
		if(y(iyp)==0)
			disp(x);
			printf("A has the eigenvalue 0, repick x and restart");
			p = x;
			return;
		endif
		ERR = infnorm(x-(y/iy));
		x = y/iy;
		if(ERR<TOL)
			printf("Iterations: %d\n",k);
			p = x;
			eigval = mu;
			return;
		endif
		k++;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	return;
endfunction

%Accelerated Power Method
function p = apm(A,x,TOL,maxIter)
	k = 1;
	muZero = 0;
	muOne = 0;
	p = 1;
	n = columns(A);
	ix = infnorm(x);
	ixp = 0;
	while(p<=n)
		if(abs(x(p))==ix)
			ixp = p;
			break;
		endif
		p++;
	endwhile
	x = x/ix;
	while(k<=maxIter)
		y = A*x;
		mu = y(ixp);
		muHat = muZero-((muOne-muZero)**2)/(mu-2*muOne+muZero);
		p = 1;
		iy = infnorm(y);
		iyp = 0;
		while(p<=n)
			if(abs(y(p))==iy)
				iyp = p;
				break;
			endif
			p++;
		endwhile
		if(y(iyp)==0)
			disp(x);
			printf("A has the eigenvalue 0, repick x and restart");
			p = x;
			return;
		endif
		ERR = infnorm(x-(y/iy));
		x = y/iy;
		if(ERR<TOL&&k>=4)
			printf("Eigenvalue:\n");
			disp(muHat);
			printf("Associated eigenvector:\n");
			disp(x);
			printf("Iterations: %d\n",k);
			p = x;
			return;
		endif
		k++;
		muZero = muOne;
		muOne = mu;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	return;
endfunction

%Inf Norm
function i = infnorm(x)
	i = max(abs(max(x)),abs(min(x)));
	return;
endfunction

%Two Norm
function t = twonorm(x)
	sum = 0;
	i = 1;
	while(i<=rows(x))
		sum= sum + x(i)**2;
		i++;
	endwhile
	t = sqrt(sum);
	return;
endfunction

%Symmetric Power Method
function p = spm(A,x,TOL,maxIter)
	if(A!=A')
		printf("Only use this for symmetric matrices\n");
		p = -1;
		return;
	endif
	k = 1;
	x = x/twonorm(x);
	while(k<=maxIter)
		y = A*x;
		mu = x'*y;
		ty = twonorm(y);
		if(ty==0)
			disp(x);
			printf("A has the eigenvalue 0, repick x and restart");
			p = x;
			return;
		endif
		ERR = twonorm(x-(y/ty));
		x = y/ty;
		if(ERR<TOL)
			printf("Eigenvalue:\n");
			disp(mu);
			printf("Associated eigenvector:\n");
			disp(x);
			printf("Iterations: %d\n",k);
			p = x;
			return;
		endif
		k++;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	return;
endfunction

%Accelerated Symmetric Power Method
%Use's Aitkens Delta-Squared Method
function p = aspm(A,x,TOL,maxIter)
	if(A!=A')
		printf("Only use this for symmetric matrices\n");
		p = -1;
		return;
	endif
	k = 1;
	x = x/twonorm(x);
	muZero = 0;
	muOne = 0;
	while(k<=maxIter)
		y = A*x;
		mu = x'*y;
		muHat = muZero-((muOne-muZero)**2)/(mu-2*muOne+muZero);
		ty = twonorm(y);
		if(ty==0)
			disp(x);
			printf("A has the eigenvalue 0, repick x and restart");
			p = x;
			return;
		endif
		ERR = twonorm(x-(y/ty));
		x = y/ty;
		if(ERR<TOL)
			printf("Eigenvalue:\n");
			disp(muOne);
			printf("Associated eigenvector:\n");
			disp(x);
			printf("Iterations: %d\n",k);
			p = x;
			return;
		endif
		k++;
		muZero = muOne;
		muOne = muHat;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	return;
endfunction

%Inverse Power Method
function p = ipm(A,x,TOL,maxIter)
	k = 1;
	p = 1;
	q = (x'*A*x)/(x'*x);
	ix = infnorm(x);
	ixp = 0;
	n = columns(A);
	while(p<=n)
		if(abs(x(p))==ix)
			ixp = p;
			break;
		endif
		p++;
	endwhile
	x = x/ix;
	I = eye(columns(A));
	while(k<=maxIter)
		y = (A-q*I)\x;
		if(y==-1)
			printf("q is an eigenvalue:\n");
			disp(q);
			p = q;
			return;
		endif
		mu = y(ixp);
		p = 1;
		iy = infnorm(y);
		iyp = 0;
		while(p<=n)
			if(abs(y(p))==iy)
				iyp = p;
				break;
			endif
			p++;
		endwhile
		ERR = infnorm(x-(y/y(iyp)));
		x = y/y(iyp);
		if(ERR<TOL)
			mu = 1/mu + q;
			printf("Eigenvalue:\n");
			disp(mu);
			printf("Associated eigenvector:\n");
			disp(x);
			printf("Iterations: %d\n",k);
			p = x;
			return;
		endif
		k++;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	return;
endfunction

%Wielandt Deflation
function [w,val] = wd(A,eigval,eigvec,x,TOL,maxIter)
	i = 1;
	n = rows(A);
	m = max(abs(eigvec));
	B = zeros(n-1);
	w = zeros(n,1);
	u = zeros(n,1);
	while(i<=rows(eigvec))
		if(abs(eigvec(i))==m)
			break;
		endif
		i++;
	endwhile
	if(i!=1)
		k = 1;
		while(k<=(i-1))
			j = 1;
			while(j<=(i-1))
				B(k,j) = A(k,j)-(eigvec(k)/eigvec(i))*A(i,j);
				j++;
			endwhile
			k++;
		endwhile
	endif
	if(i!=1&&i!=n)
		k = i;
		while(k<=(n-1))
			j = 1;
			while(j<=(i-1))
				B(k,j) = A(k+1,j)-(eigvec(k+1)/eigvec(i))*A(i,j);
				B(j,k) = A(j,k+1)-(eigvec(j)/eigvec(i))*A(i,k+1);
				j++;
			endwhile
			k++;
		endwhile
	endif
	if(i!=n)
		k = i;
		while(k<=(n-1))
			j = i;
			while(j<=(n-1))
				B(k,j) = A(k+1,j+1)-(eigvec(k+1)/eigvec(i))*A(i,j+1);
				j++;
			endwhile
			k++;
		endwhile
	endif
	[wPrime,mu] = pm(B,x,TOL,maxIter);
	if(wPrime==-1&&mu==-1)
		printf("Failed");
		w = -1;
		return;
	endif
	if(i!=1)
		k = 1;
		while(k<=(i-1))
			w(k) = wPrime(k);
			k++;
		endwhile
	endif
	w(i)=0;
	if(i!=n)
		k = i+1;
		while(k<=n)
			w(k) = wPrime(k-1);
			k++;
		endwhile
	endif
	k = 1;
	while(k<=n)
		j = 1;
		sum = 0;
		while(j<=n)
			sum = sum + A(i,j)*w(j);
			j++;
		endwhile
		u(k) = (mu-eigval)*w(k)+(eigvec(k)/eigvec(i))*sum;
		k++;
	endwhile
	w = u;
	val = mu;
	return;
endfunction

%Inverse Power Method for Wielandt Deflation
function [p,v] = ipmw(A,x,TOL,maxIter)
	k = 1;
	p = 1;
	q = (x'*A*x)/(x'*x);
	ix = infnorm(x);
	ixp = 0;
	n = columns(A);
	while(p<=n)
		if(abs(x(p))==ix)
			ixp = p;
			break;
		endif
		p++;
	endwhile
	x = x/ix;
	I = eye(columns(A));
	while(k<=maxIter)
		y = (A-q*I)\x;
		val = errno();
		if(val==9)
			printf("q is an eigenvalue:\n");
			disp(q);
			p = q;
			v = x;
			return;
		endif
		mu = y(ixp);
		p = 1;
		iy = infnorm(y);
		iyp = 0;
		while(p<=n)
			if(abs(y(p))==iy)
				iyp = p;
				break;
			endif
			p++;
		endwhile
		ERR = infnorm(x-(y/y(iyp)));
		x = y/y(iyp);
		if(ERR<TOL)
			mu = 1/mu + q;
			p = x;
			v = mu;
			return;
		endif
		k++;
	endwhile
	printf("Max number of iterations exceeded\n");
	p = -1;
	v = -1;
	return;
endfunction

%Householder's Method
function h = hh(A)
	n = rows(A);
	k = 1;
	while(k<=(n-2))
		v = zeros(1,n);
		u = zeros(1,n);
		z = zeros(1,n);
		j = k+1;
		q = 0;
		while(j<=n)
			q = q + A(j,k)**2;
			j++;
		endwhile
		if(A(k+1,k)==0)
			alpha = -sqrt(q);
		else
			alpha = -(sqrt(q)*A(k+1,k))/(abs(A(k+1,k)));
		endif
		RSQ = alpha**2-alpha*A(k+1,k);
		v(k) = 0;
		v(k+1) = A(k+1,k)-alpha;
		j = k+2;
		while(j<=n)
			v(j) = A(j,k);
			j++;
		endwhile
		j = k;
		while(j<=n)
			i = k+1;
			sum = 0;
			while(i<=n)
				sum = sum + A(j,i)*v(i);
				i++;
			endwhile
			u(j) = sum*(1/RSQ);
			j++;
		endwhile
		i = k+1;
		PROD = 0;
		while(i<=n)
			PROD = PROD + u(i)*v(i);
			i++;
		endwhile
		j = k;
		while(j<=n)
			z(j) = u(j)-(PROD/(2*RSQ))*v(j);
			j++;
		endwhile
		l = k+1;
		while(l<=n-1)
			j = l+1;
			while(j<=n)
				A(j,l) = A(j,l)-v(l)*z(j)-v(j)*z(l);
				A(l,j) = A(j,l);
				j++;
			endwhile
			A(l,l) = A(l,l)-2*v(l)*z(l);
			l++;
		endwhile
		A(n,n) = A(n,n)-2*v(n)*z(n);
		j = k+2;
		while(j<=n)
			A(k,j) = 0;
			A(j,k) = 0;
			j++;
		endwhile
		A(k+1,k) = A(k+1,k)-v(k+1)*z(k);
		A(k,k+1) = A(k+1,k);
		k++;
	endwhile
	h = A;
	return;
endfunction

%Upper Hessenberg Variation of Householder's Method
function h = upperh(A)
	n = rows(A);
	k = 1;
	while(k<=(n-2))
		v = zeros(1,n);
		u = zeros(1,n);
		z = zeros(1,n);
		y = zeros(1,n);
		j = k+1;
		q = 0;
		while(j<=n)
			q = q + A(j,k)**2;
			j++;
		endwhile
		if(A(k+1,k)==0)
			alpha = -sqrt(q);
		else
			alpha = -(sqrt(q)*A(k+1,k))/(abs(A(k+1,k)));
		endif
		RSQ = alpha**2-alpha*A(k+1,k);
		v(k) = 0;
		v(k+1) = A(k+1,k)-alpha;
		j = k+2;
		while(j<=n)
			v(j) = A(j,k);
			j++;
		endwhile

		%Step 6
		j = 1;
		while(j<=n)
			i = k+1;
			sum = 0;
			sumTwo = 0;
			while(i<=n)
				sum = sum + A(j,i)*v(i);
				sumTwo = sumTwo + A(i,j)*v(i);
				i++;
			endwhile
			u(j) = sum*(1/RSQ);
			y(j) = sumTwo*(1/RSQ);
			j++;
		endwhile
		
		%Step 7
		i = k+1;
		PROD = 0;
		while(i<=n)
			PROD = PROD + u(i)*v(i);
			i++;
		endwhile
		
		%Step 8
		j = 1;
		while(j<=n)
			z(j) = u(j)-(PROD/RSQ)*v(j);
			j++;
		endwhile
		
		%Step 9
		l = k+1;
		while(l<=n)
			%Step 10
			j = 1;
			while(j<=k)
				A(j,l) = A(j,l)-v(l)*z(j);
				A(l,j) = A(l,j)-y(j)*v(l);
				j++;
			endwhile
			
			%Step 11
			j = k+1;
			while(j<=n)
				A(j,l) = A(j,l)-z(j)*v(l)-y(l)*v(j);
				j++;
			endwhile
			l++;
		endwhile
		k++;
	endwhile
	h = A;
	return;
endfunction

#QR Algorithm
function qq = qr(A,B,TOL, maxIter)
	k = 1;
	n = columns(A);
	shift = 0;
	while(k<=maxIter)
		%Step 3
		if(abs(B(n))<=TOL)
			lambda = A(n) + shift;
			disp(lambda);
			n--;
		endif
		%Step 4
		if(abs(B(2))<=TOL)
			lambda = A(1) + shift;
			disp(lambda);
			n--;
			A(1) = A(2);
			j = 2;
			while(j<=n)
				A(j) = A(j+1);
				B(j) = B(j+1);
				j++;
			endwhile
		endif
		%Step 5
		if(n==0) 
			return;
		endif
		%Step 6
		if(n==1)
			lambda = A(1) + shift;
			disp(lambda);
			return;
		endif
		%Step 7
		j = 3;
		while(j<=(n-1))
			if(abs(B(j))<=TOL)
				printf("Split matrix at %d\n",j);
				return;
			endif
			j++;
		endwhile
		%Step 8
		b = -(A(n-1) + A(n));
		c = A(n)*A(n-1)-B(n)**2;
		d = sqrt(b**2-4*c);
		%Step 9
		if(b>0)
			muOne = (-2*c)/(b+d);
			muTwo = -(b+d)/2;
		else
			muOne = (d-b)/2;
			muTwo = (2*c)/(d-b);
		endif
		%Step 10
		if(n==2)
			lambda = muOne + shift;
			disp(lambda);
			lambda = muTwo + shift;
			disp(lambda);
			return;
		endif
		%Step 11
		if(abs(muOne-A(n)) <= abs(muTwo-A(n)))
			delta = muOne;
		else
			delta = muTwo;
		endif
		%Step 12
		shift = shift + delta;
		j = 1;
		d = zeros(1,n);
		%Step 13
		while(j<=n)
			d(j) = A(j) - delta;
			j++;
		endwhile
		x = zeros(1,n);
		y = zeros(1,n);
		delta = zeros(1,n);
		c = zeros(1,n);
		z = zeros(1,n);
		%Step 14
		x(1) = d(1);
		y(1) = B(2);
		q = zeros(1,n);
		r = zeros(1,n);
		%Step 15
		j = 2;
		while(j<=n)
			z(j-1) = sqrt(x(j-1)**2 + B(j)**2);
			c(j) = x(j-1)/z(j-1);
			delta(j) = B(j)/z(j-1);
			q(j-1) = c(j)*y(j-1)+delta(j)*d(j);
			x(j) = -delta(j)*y(j-1)+c(j)*d(j);
			if(j!=n)
				r(j-1) = delta(j)*B(j+1);
				y(j) = c(j)*B(j+1);
			endif
			j++;
		endwhile
		%Step 16
		z(n) = x(n);
		A(1) = delta(2)*q(1)+c(2)*z(1);
		B(2) = delta(2)*z(2);
		%Step 17
		j = 2;
		while(j<=(n-1))
			A(j) = delta(j+1)*q(j)+c(j)*c(j+1)*z(j);
			B(j+1) = delta(j+1)*z(j+1);
			j++;
		endwhile
		%Step 18
		A(n) = c(n)*z(n);
		k++;
	endwhile
	printf("Max iterations exceeded");
	qq = -1;
	return;
endfunction