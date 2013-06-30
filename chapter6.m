#Gaussian Elimination
function g = gauss(mat)
	n = rows(mat);
	g = zeros(1,n);
	for i = 1:n
		j = i;
		p = 0;
		while(j<=n)
			if(mat(j,i)!=0)
				p = j;
				break;
			endif
			j++;
		endwhile
		if(p==0)
			disp("No unique solution exists");
			g=-1;
			return;
		elseif(p!=i)
			printf("Row interchange: %3d  with %3d\n",i,p);
			temp = mat(p,:);
			mat(p,:)=mat(i,:);
			mat(i,:)=temp;
		endif
		for k = (i+1):n
			mki = mat(k,i)/mat(i,i);
			mat(k,:) = mat(k,:) - mki*mat(i,:);
		endfor
	endfor
	if(mat(n,n)==0)
		disp("No unique solution exists");
		g=-1;
		return;
	endif
	g(n) = mat(n,n+1)/mat(n,n);
	i = n-1;
	while(i>=1)
		sum = 0;
		for j = i+1:n
			sum = sum + mat(i,j)*g(j);
		endfor
		g(i) = (mat(i,n+1)-sum)/mat(i,i);
		i--;
	endwhile
	return;
endfunction

#Gauss-Jordan Elimination
function g = gj(mat)
	n = rows(mat);
	g = zeros(1,n);
	for i = 1:n
		j = i;
		p = 0;
		while(j<=n)
			if(mat(j,i)!=0)
				p = j;
				break;
			endif
			j++;
		endwhile
		if(p==0)
			disp("No unique solution exists");
			g=-1;
			return;
		elseif(p!=i)
			printf("Row interchange: %3d  with %3d\n",i,p);
			temp = mat(p,:);
			mat(p,:)=mat(i,:);
			mat(i,:)=temp;
		endif
		for k = (i+1):n
			mki = mat(k,i)/mat(i,i);
			mat(k,:) = mat(k,:) - mki*mat(i,:);
		endfor
	endfor
	if(mat(n,n)==0)
		disp("No unique solution exists");
		g=-1;
		return;
	endif
	i = n;
	while(i>=1)
		j = i-1;
		while(j>=1)
			mult = mat(j,i)/mat(i,i);
			mat(j,i) = 0;
			mat(j,n+1) =  mat(j,n+1) - mat(i,n+1)*mult;
			j--;
		endwhile
		i--;
	endwhile
	for j = 1:n
		mat(j,n+1) = mat(j,n+1)/mat(j,j);
		mat(j,j) = 1;
	endfor
	g = mat;
	return;
endfunction

#Gauss Elimination with Partial Pivoting
function g = gpp(mat)
	n = rows(mat);
	g = zeros(1,n);
	for i = 1:n
		j = i;
		p = 0;
		while(j<=n)
			if(p==0||abs(mat(j,i)) > abs(mat(p,i)))
				p = j;
			endif
			j++;
		endwhile
		if(p==0)
			disp("No unique solution exists");
			g=-1;
			return;
		elseif(p!=i)
			printf("Row interchange: %3d  with %3d\n",i,p);
			temp = mat(p,:);
			mat(p,:)=mat(i,:);
			mat(i,:)=temp;
		endif
		for k = (i+1):n
			mki = mat(k,i)/mat(i,i);
			mat(k,:) = mat(k,:) - mki*mat(i,:);
		endfor
	endfor
	if(mat(n,n)==0)
		disp("No unique solution exists");
		g=-1;
		return;
	endif
	g(n) = mat(n,n+1)/mat(n,n);
	i = n-1;
	while(i>=1)
		sum = 0;
		for j = i+1:n
			sum = sum + mat(i,j)*g(j);
		endfor
		g(i) = (mat(i,n+1)-sum)/mat(i,i);
		i--;
	endwhile
	return;
endfunction

#Gauss Elimination with Scaled Partial Pivoting
function g = gspp(mat)
	n = rows(mat);
	g = zeros(1,n);
	s = zeros(1,n);
	for i = 1:n
		for j = 1:n
			if(abs(mat(i,j)) > s(i))
				s(i) = abs(mat(i,j));
			endif
		endfor
		if(s(i)==0)
			disp("No unique solution exists");
			break;
		endif
	endfor
	for i = 1:n
		j = i;
		p = 0;
		while(j<=n)
			if(p==0||(abs(mat(j,i))/s(j) > abs(mat(p,i))/s(p)))
				p = j;
			endif
			j++;
		endwhile
		if(p==0)
			disp("No unique solution exists");
			g=-1;
			return;
		elseif(p!=i)
			printf("Row interchange: %3d  with %3d\n",i,p);
			temp = mat(p,:);
			mat(p,:)=mat(i,:);
			mat(i,:)=temp;
		endif
		for k = (i+1):n
			mki = mat(k,i)/mat(i,i);
			mat(k,:) = mat(k,:) - mki*mat(i,:);
		endfor
	endfor
	if(mat(n,n)==0)
		disp("No unique solution exists");
		g=-1;
		return;
	endif
	g(n) = mat(n,n+1)/mat(n,n);
	i = n-1;
	while(i>=1)
		sum = 0;
		for j = i+1:n
			sum = sum + mat(i,j)*g(j);
		endfor
		g(i) = (mat(i,n+1)-sum)/mat(i,i);
		i--;
	endwhile
	return;
endfunction

#LU Factorization, 1s on Diagonal of Lower Matrix
function l = lu(mat)
	n = rows(mat);
	lower = zeros(n);
	upper = zeros(n);
	if(n!=columns(mat))
		disp("Matrix must be square");
		return;
	elseif(mat(1,1)==0)
		disp("Factorization impossible");
		return;
	else
		lower(1,1) = 1;
		upper(1,1) = mat(1,1);
	endif
	for j = 2:n
		upper(1,j) = mat(1,j)/lower(1,1);
		lower(j,1) = mat(j,1)/upper(1,1);
	endfor
	for i = 2:(n-1)
		sum  = mat(i,i);
		for k = 1:i-1
			sum = sum - lower(i,k)*upper(k,i);
		endfor
		if(sum == 0)
			disp("Factorization impossible");
			return;
		else
			lower(i,i) = 1;
			upper(i,i) = sum;
		endif
		for j = (i+1):n
			sumup = mat(i,j);
			sumlow = mat(j,i);
			for k = 1:(i-1)
				sumup = sumup - lower(i,k)*upper(k,j);
				sumlow = sumlow - lower(j,k)*upper(k,i);
			endfor
			upper(i,j) = (1/lower(i,i))*sumup;
			lower(j,i) = (1/upper(i,i))*sumlow;
		endfor
	endfor
	sumlast = mat(n,n);
	for k = 1:(n-1)
		sumlast = sumlast - lower(n,k)*upper(k,n);
	endfor
	lower(n,n) = 1;
	upper(n,n) = sumlast;
	printf("Lower:\n");
	disp(lower);
	printf("Upper:\n");
	disp(upper);
endfunction

#LDLt Factorization
function l = ldlt(mat)
	n = rows(mat);
	lower = zeros(n);
	diagonal = zeros(n);
	v = zeros(1,n);
	if(n!=columns(mat))
		disp("Matrix must be square");
		return;
	endif
	for i = 1:n
		sum = mat(i,i);
		lower(i,i) = 1;
		for j = 1:i-1
			v(j) = lower(i,j)*diagonal(j,j);
			sum = sum - lower(i,j)*v(j);
		endfor
		diagonal(i,i) = sum;
		for j = (i+1):n
			sumlv = mat(j,i);
			for k = 1:(i-1)
				sumlv = sumlv - lower(j,k)*v(k);
			endfor
			lower(j,i) = sumlv/diagonal(i,i);
		endfor
	endfor
	printf("Lower:\n");
	disp(lower);
	printf("Diagonal:\n");
	disp(diagonal);
	printf("Lower Transpose:\n");
	disp(lower.');
endfunction

#Cholesky Factorization
function l = cho(mat)
	n = rows(mat);
	lower = zeros(n);
	if(n!=columns(mat))
		disp("Matrix must be square");
		return;
	endif
	lower(1,1) = sqrt(mat(1,1));
	for j = 2:n
		lower(j,1) = mat(j,1)/lower(1,1);
	endfor
	for i = 2:(n-1)
		sum = mat(i,i);
		for k = 1:(i-1)
			sum = sum - lower(i,k)*lower(i,k);
		endfor
		lower(i,i) = sqrt(sum);
		for j = (i+1):n
			sumj = mat(j,i);
			for k = 1:(i-1)
				sumj = sumj - lower(j,k)*lower(i,k);
			endfor
			lower(j,i) = (1/lower(i,i))*sumj;
		endfor
	endfor
	sumn = mat(n,n);
	for k = 1:(n-1)
		sumn = sumn - lower(n,k)*lower(n,k);
	endfor
	lower(n,n) = sqrt(sumn);
	printf("Lower:\n");
	disp(lower);
	printf("Lower Transpose:\n");
	disp(lower.');
endfunction

#Crout Factorization for Tridiagonal Linear Systems
function lo = crout(mat)
	n = rows(mat);
	l = zeros(n);
	z = zeros(1,n);
	x = zeros(1,n);
	u = zeros(n);
	l(1,1) = mat(1,1);
	u(1,2) = mat(1,2)/l(1,1);
	z(1) = mat(1,n+1)/l(1,1);
	for i = 2:n-1
		l(i,i-1) = mat(i,i-1);
		l(i,i) = mat(i,i)-l(i,i-1)*u(i-1,i);
		u(i,i+1) = mat(i,i+1)/l(i,i);
		z(i) = (mat(i,n+1) - l(i,i-1)*z(i-1))/l(i,i);
	endfor
	l(n,n-1) = mat(n,n-1);
	l(n,n) = mat(n,n) - l(n,n-1)*u(n-1,n);
	z(n) = (mat(n,n+1)-l(n,n-1)*z(n-1))/l(n,n);
	x(n) = z(n);
	i = n-1;
	while(i>=1)
		x(i) = z(i) - u(i,i+1)*x(i+1);
		i--;
	endwhile
	printf("Solution:");
	disp(x);
endfunction