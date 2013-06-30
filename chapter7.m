#Jacobi Iteration
function j = jacobi(A,b,x,TOL)
	newx = x;
	it = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sum = 0;
			for k = 1:columns(x)
				if(k!=i)
					sum = sum + (-A(i,k)*x(k));
				endif
			endfor
			sum = sum + b(i);
			sum = sum/A(i,i);
			newx(i) = sum;
		endfor
		it++;
	until (infnorm(newx-x)/infnorm(newx) < TOL)
	printf("Iterations: %d\n",it);
	j = newx;
	return;
endfunction

#Jacobi Iteration with Iteration Bound
function j = jacobiit(A,b,x,it)
	newx = x;
	iteration = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sum = 0;
			for k = 1:columns(x)
				if(k!=i)
					sum = sum + (-A(i,k)*x(k));
				endif
			endfor
			sum = sum + b(i);
			sum = sum/A(i,i);
			newx(i) = sum;
		endfor
		iteration++;
		disp(newx);
	until (iteration>=it)
	j = newx;
	return;
endfunction

function i = infnorm(x)
	i = max(abs(max(x)),abs(min(x)));
	return;
endfunction

#Gauss-Seidel Iteration
function g = gsi(A,b,x,TOL)
	newx = x;
	it = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sumOne = 0;
			sumTwo = 0;
			for k = 1:(i-1)
				sumOne = sumOne + (A(i,k)*newx(k));
			endfor
			for k = (i+1):columns(x)
				sumTwo = sumTwo + (A(i,k)*x(k));
			endfor
			sum = -sumOne-sumTwo + b(i);
			sum = sum/A(i,i);
			newx(i) = sum;
		endfor
		it++;
	until (infnorm(newx-x)/infnorm(newx) < TOL)
	printf("Iterations: %d\n",it);
	g = newx;
	return;
endfunction

#Gauss-Seidel Iteration with Iteration Bound
function g = gsib(A,b,x,TOL, iter)
	newx = x;
	it = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sumOne = 0;
			sumTwo = 0;
			for k = 1:(i-1)
				sumOne = sumOne + (A(i,k)*newx(k));
			endfor
			for k = (i+1):columns(x)
				sumTwo = sumTwo + (A(i,k)*x(k));
			endfor
			sum = -sumOne-sumTwo + b(i);
			sum = sum/A(i,i);
			printf("sum: %d\n",sum);
			newx(i) = sum;
		endfor
		it++;
	until ((infnorm(newx-x)/infnorm(newx) < TOL)||(it>=iter))
	if(it>=iter)
		printf("Failed");
		return;
	endif
	printf("Iterations: %d\n",it);
	g = newx;
	return;
endfunction

#Successive Over-Relaxation
function g = sor(A,b,x,omega,TOL)
	newx = x;
	it = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sumOne = 0;
			sumTwo = 0;
			for k = 1:(i-1)
				sumOne = sumOne + (A(i,k)*newx(k));
			endfor
			for k = (i+1):columns(x)
				sumTwo = sumTwo + (A(i,k)*x(k));
			endfor
			sum = -sumOne-sumTwo + b(i);
			sum = omega*sum/A(i,i);
			newx(i) = (1-omega)*x(i)+sum;
		endfor
		it++;
	until (infnorm(newx-x)/infnorm(newx) < TOL)
	printf("Iterations: %d\n",it);
	g = newx;
	return;
endfunction

#Successive Over-Relaxation with Iteration Bound
function g = sorb(A,b,x,omega,TOL,bound)
	newx = x;
	it = 0;
	do
		x = newx;
		for i = 1:columns(x)
			sumOne = 0;
			sumTwo = 0;
			for k = 1:(i-1)
				sumOne = sumOne + (A(i,k)*newx(k));
			endfor
			for k = (i+1):columns(x)
				sumTwo = sumTwo + (A(i,k)*x(k));
			endfor
			sum = -sumOne-sumTwo + b(i);
			sum = omega*sum/A(i,i);
			newx(i) = (1-omega)*x(i)+sum;
		endfor
		it++;
		disp(newx);
	until ((infnorm(newx-x)/infnorm(newx)) < TOL||(it>=bound))
	printf("Iterations: %d\n",it);
	g = newx;
	return;
endfunction

#Inf Norm of a Matrix
function a = infnorm(mat)
	largest = 0;
	for i = 1:rows(mat)
		sum = 0;
		for j = 1:columns(mat)
			sum+=abs(mat(i,j));
		endfor
		if(sum > largest)
			largest = sum;
		endif
	endfor
	a = largest;
	return;
endfunction

#Condition of a Matrix
function c = matcond(mat)
	c = infnorm(mat)*infnorm(inv(mat));
	return;
endfunction

#Negative Square Diagonal Matrix
function d = dneg(A)
	d = zeros(rows(A));
	for i = 1:rows(A)
		d(i,i) = 1/sqrt(A(i,i));
	endfor
	return;
endfunction

#Conjugate Gradient Method
function c = pcg(A,b,C,x,maxIter,TOL)
	n = rows(A);
	r = b-A*x;
	w = C*r;
	v = C*w;
	j = 1;
	alpha = 0;
	while(j<=n)
		alpha = alpha + w(j)**2;
		j++;
	endwhile
	k = 1;
	while(k<=maxIter)
		if(infnorm(v)<TOL)
			printf("Residual vector r:\n");
			disp(r);
			printf("Number of iterations: %d\n", k);
			c = x;
			return;
		endif
		u = A*v;
		j = 1;
		sum = 0;
		while(j<=n)
			sum = sum + v(j)*u(j);
			j++;
		endwhile
		t = alpha/sum;
		x = x+t*v;
		r = r-t*u;
		w = C*r;
		j = 1;
		sum = 0;
		while(j<=n)
			sum = sum + w(j)**2;
			j++;
		endwhile
		if((abs(sum) < TOL)&&(infnorm(r)<TOL))
			printf("Residual vector r:\n");
			disp(r);
			printf("Number of iterations: %d\n", k);
			c = x;
			return;
		endif
		s = sum/alpha;
		v = C*w+s*v;
		alpha = sum;
		k++;
	endwhile
	printf("Procedure unsuccessful, max iterations exceeded\n");
	c = x;
	return;
endfunction

function s = sevensix(i)
	if(i==0)
		s = zeros(25);
		for i = 1:25
			for j = 1:25
				if(i==j)
					s(i,j)=4;
				elseif((j==i+1)&&(i<=24))
					s(i,j) = -1;
				elseif((j==i-1)&&(i>=2))
					s(i,j) = -1;
				elseif((j==i+5)&&(i<=20))
					s(i,j) = -1;
				elseif((j==i-5)&&(i>=6))
					s(i,j) = -1;
				endif
			endfor
		endfor
	elseif(i==1)
		s = zeros(40);
		for i = 1:40
			for j = 1:40
				if(i==j)
					s(i,j)=2*i;
				elseif((j==i+1)&&(i<=39))
					s(i,j) = -1;
				elseif((j==i-1)&&(i>=2))
					s(i,j) = -1;
				endif
			endfor
		endfor
	endif
	return;
endfunction