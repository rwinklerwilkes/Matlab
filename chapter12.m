%Poisson Equation Finite-Difference
%Variable key
%a - Lowest x-value
%b - Highest x-value
%c - Lowest y-value
%d - Highest y-value
%m - Number of y partitions, >=3
%n - Number of x partitions, >=3
%f - Function they satisfy
%gB - Boundary condition for y=c
%gT - Boundary condition for y=d
%gL - Boundary condition for x=a
%gR - Boundary condition for x=b
function p = pefd(a,b,c,d,m,n,f,gB,gT,gL,gR,TOL,maxIter)
	%Step 1
	h = (b-a)/n;
	k = (d-c)/m;
	x = zeros(1,n-1);
	y = zeros(1,m-1);
	%Step 4
	w = zeros(n,m);
	%Step 2
	i=1;
	while(i<=(n-1))
		x(i) = a+i*h;
		i++;
	endwhile
	%Step 3
	i = 1;
	while(i<=(m-1))
		y(i) = c+i*k;
		i++;
	endwhile
	%Step 5
	lambda = (h**2)/(k**2);
	mu = 2*(1+lambda);
	l = 1;
	%Step 6
	while(l<=maxIter)
		%Step 7
		z = ((-h**2)*f(x(1),y(m-1)) +gL(y(m-1)) + lambda*gT(x(1)) + lambda*w(1,m-2) + w(2,m-1))/mu;
		NORM = abs(z-w(1,m-1));
		w(1,m-1) = z;
		%Step 8
		i = 2;
		while(i<=(n-2))
			z = ((-h**2)*f(x(i),y(m-1)) + lambda*gT(x(i)) + w(i-1,m-1) + w(i+1,m-1) + lambda*w(i,m-2))/mu;
			if(abs(w(i,m-1) - z) > NORM)
				NORM = abs(w(i,m-1)-z);
			endif
			w(i,m-1) = z;
			i++;
		endwhile
		%Step 9
		z = ((-h**2)*f(x(n-1),y(m-1)) + gR(y(m-1)) + lambda*gT(x(n-1)) + w(n-2,m-1) + lambda*w(n-1,m-2))/mu;
		if(abs(w(n-1,m-1)-z) > NORM)
			NORM = abs(w(n-1,m-1)-z);
		endif
		w(n-1,m-1) = z;
		%Step 10
		j = m-2;
		while(j>=2)
			%Step 11
			z = ((-h**2)*f(x(1),y(j)) + gL(y(j)) + lambda*w(1,j+1)+lambda*w(1,j-1)+w(2,j))/mu;
			if(abs(w(1,j)-z) > NORM)
				NORM = abs(w(1,j)-z);
			endif
			w(1,j) = z;
			%Step 12
			i = 2;
			while(i<=n-2)
				z = ((-h**2)*f(x(i),y(j)) + w(i-1,j) + lambda*w(i,j+1)+w(i+1,j)+lambda*w(i,j-1))/mu;
				if(abs(w(i,j)-z) > NORM)
					NORM = abs(w(i,j)-z);
				endif
				w(i,j) = z;
				i++;
			endwhile
			%Step 13
			z= (-(h**2)*f(x(n-1),y(j)) + gR(y(j)) + w(n-2,j) + lambda*w(n-1,j+1) + lambda*w(n-1,j-1))/mu;
			if(abs(w(n-1,j)-z) > NORM)
				NORM = abs(w(n-1,j)-z);
			endif
			w(n-1,j) = z;
			j--;
		endwhile
		%Step 14
		z= (-(h**2)*f(x(1),y(1)) + gL(y(1))+lambda*gB(x(1)) + lambda*w(1,2) + w(2,1))/mu;
		if(abs(w(1,1)-z) > NORM)
			NORM = abs(w(1,1)-z);
		endif
		w(1,1) = z;
		%Step 15
		i = 2;
		while(i<=(n-2))
			z= (-(h**2)*f(x(i),y(1)) + lambda*gB(x(i))+w(i-1,1) + lambda*w(i,2) + w(i+1,1))/mu;
			if(abs(w(i,1)-z) > NORM)
				NORM = abs(w(i,1)-z);
			endif
			w(i,1) = z;
			i++;
		endwhile
		%Step 16
		z = (-(h**2)*f(x(n-1),y(1)) + gR(y(1)) + lambda*gB(x(n-1)) + w(n-2,1) + lambda*w(n-1,2))/mu;
		if(abs(w(n-1,1) - z) > NORM)
			NORM = abs(w(n-1,1) - z);
		endif
		w(n-1,1) = z;
		%Step 17
		if(NORM<=TOL)
			i = 1;
			while(i<=(n-1))
				j = 1;
				while(j<=(m-1))
					printf("x(%d): %g, y(%d): %g, w(%d,%d): %g\n", i,x(i),j,y(j),i,j,w(i,j));
					j++;
				endwhile
				i++;
			endwhile
			p = 0;
			return;
		endif
		l++;
	endwhile
	printf("Max iterations exceeded");
	p = -1;
	return;
endfunction

%g - Dirichlet boundary condition
%g1,g2 - Neumann boundary conditions
%verts - Vertices of the triangulation
%nodes - each node corresponds to the row of verts which contains the vertex for that node
%triVerts - each row has 3 values - the three node numbers of the triangle
%K - Interior triangles
%N - Triangles with at least one Neumann edge
%M - All other triangles (total triangles - K+N+M)
%n - Nodes either in D or on the Neumann edge
%m - Dirichlet nodes
function f = fe(p,q,r,f,g,g1,g2,verts,nodes,triVerts,K,N,M,n,m)
	gamma = zeros(1,m);
	l = n+1;
	while(l<=m)
		gamma(l) = g(verts(nodes(l),1),verts(nodes(l),2));
		l++;
	endwhile
	alpha = zeros(n,n);
	beta = zeros(1,n);
	delta = zeros(1,M);
	a = zeros(M,3);
	b = zeros(M,3);
	c = zeros(M,3);
	z = zeros(M,3,3);
	J = zeros(M,3,3);
	I = zeros(M,3);
	H = zeros(M,3);
	i = 1;
	while(i<=M)
		x1 = verts(nodes(triVerts(i,1)),1);
		x2 = verts(nodes(triVerts(i,2)),1);
		x3 = verts(nodes(triVerts(i,3)),1);
		y1 = verts(nodes(triVerts(i,1)),2);
		y2 = verts(nodes(triVerts(i,2)),2);
		y3 = verts(nodes(triVerts(i,3)),2);
		delta(i) = det([1,x1,y1;1,x2,y2;1,x3,y3]);
		a(i,1) = (x2*y3-y2*x3)/delta(i);
		a(i,2) = (x3*y1-y3*x1)/delta(i);
		a(i,3) = (x1*y2-y1*x2)/delta(i);
		b(i,1) = (y2-y3)/delta(i);
		b(i,2) = (y3-y1)/delta(i);
		b(i,3) = (y1-y2)/delta(i);
		c(i,1) = (x3-x2)/delta(i);
		c(i,2) = (x1-x3)/delta(i);
		c(i,3) = (x2-x1)/delta(i);
		i++;
	endwhile
	i = 1;
	while(i<=M)
		j = 1;
		while(j<=3)
			k = 1;
			while(k<=j)
				nj = @(x,y)(a(i,j) + b(i,j)*x + c(i,j)*y);
				nk = @(x,y)(a(i,k) + b(i,k)*x + c(i,k)*y);
				%Calculate the integrals over the triangles
				zIntegral = @(x,y)(r(x,y)*nj(x,y)*nk(x,y));
				%Approximate the integral using quadrature rule for triangles
				%node 1
				x1 = verts(nodes(triVerts(i,1)),1);
				y1 = verts(nodes(triVerts(i,1)),2);
				%node 2
				x2 = verts(nodes(triVerts(i,2)),1);
				y2 = verts(nodes(triVerts(i,2)),2);
				%node 3
				x3 = verts(nodes(triVerts(i,3)),1);
				y3 = verts(nodes(triVerts(i,3)),2);
				area = (1/6)*abs(det([x1,x2,x3;y1,y2,y3;1,1,1]));
				%mid1 = midpoint between node 1 and node 2
				midx1 = (x1+x2)/2;
				midy1 = (y1+y2)/2;
				%mid1 = midpoint between node 2 and node 3
				midx2 = (x2+x3)/2;
				midy2 = (y2+y3)/2;
				%mid1 = midpoint between node 1 and node 3
				midx3 = (x1+x3)/2;
				midy3 = (y1+y3)/2;
				firstTerm = b(i,j)*b(i,k)*area*(p(x1,y1)+p(x2,y2)+p(x3,y3));
				secondTerm = c(i,j)*c(i,k)*area*(q(x1,y1)+q(x2,y2)+q(x3,y3));
				thirdTerm = -1*area*(zIntegral(x1,y1)+zIntegral(x2,y2)+zIntegral(x3,y3));
				z(i,j,k) = firstTerm+secondTerm+thirdTerm;
				hInt = @(x,y)(f(x,y)*nj(x,y));
				H(i,j) = -1*area*(hInt(x1,y1)+hInt(x2,y2)+hInt(x3,y3));
				k++;
			endwhile
			j++;
		endwhile
		i++;
	endwhile
	%Approximate the line integrals for the Neumann nodes
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Unfinished, probably doesn't work
	i = K+1;
	while(i<N)
		j = 1;
		while(j<=3)
			k = 1;
			while(k<=j)
				node1X = verts(nodes(i),1);
				node1Y = verts(nodes(i),2);
				node2X = verts(nodes(i+1),1);
				node2Y = verts(nodes(i+1),2);
				length = sqrt((node2X-node1X)**2+(node2Y-node1Y)**2);
				nj = @(x,y)(a(i,j) + b(i,j)*x + c(i,j)*y);
				nk = @(x,y)(a(i,k) + b(i,k)*x + c(i,k)*y);
				firstInt = @(x,y)(g1(x,y)*nj(x,y)*nk(x,y));
				secondInt = @(x,y)(g2(x,y)*nj(x,y));
				J(i,j,k) = (length/4)*(firstInt(node1X,node1Y)+firstInt(node2X,node2Y));
				I(i,j) = (length/4)*(secondInt(node1X,node1Y)+secondInt(node2X,node2Y));
				k++;
			endwhile
			j++;
		endwhile
		i++;
	endwhile
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Step 7
	i = 1;
	while(i<=M)
		k = 1;
		while(k<=3)
			%Step 8
			triNodeIKX = verts(nodes(triVerts(i,k)),1);
			triNodeIKY = verts(nodes(triVerts(i,k)),2);
			l = 1;
			while(l<=M)
				if((verts(nodes(l),1)==triNodeIKX)&&(verts(nodes(l),2)==triNodeIKY))
					break;
				endif
				l++;
			endwhile
			%Step 9
			if(k>1)
				j = 1;
				while(j<=k-1)
					%Step 10
					triNodeIKX = verts(nodes(triVerts(i,j)),1);
					triNodeIKY = verts(nodes(triVerts(i,j)),2);
					t = 1;
					while(t<=M)
						if((verts(nodes(t),1)==triNodeIKX)&&(verts(nodes(t),2)==triNodeIKY))
							break;
						endif
						t++;
					endwhile
					%Step 11
					if(l<=n)
						if(t<=n)
							alpha(l,t) = alpha(l,t)+z(i,k,j);
							alpha(t,l) = alpha(t,l)+z(i,k,j);
						else
							beta(l) = beta(l) - gamma(t)*z(i,k,j);
						endif
					else
						if(t<=n)
							beta(t) = beta(t) - gamma(l)*z(i,k,j);
						endif
					endif
					j++;
				endwhile
			endif
			if(l<=n)
				alpha(l,l) = alpha(l,l) + z(i,k,k);
				beta(l) = beta(l) + H(i,k);
			endif
			k++;
		endwhile
		i++;
	endwhile
	%Step 13
	i = K+1;
	while(i<=N)
		%Step 14
		k = 1;
		while(k<=3)
			%Step 15
			triNodeIKX = verts(nodes(triVerts(i,k)),1);
			triNodeIKY = verts(nodes(triVerts(i,k)),2);
			l = 1;
			while(l<=M)
				if((verts(nodes(l),1)==triNodeIKX)&&(verts(nodes(l),2)==triNodeIKY))
					break;
				endif
				l++;
			endwhile
			%Step 16
			if(k>1)
				j = 1;
				while(j<=k-1)
					%Step 17
					triNodeIKX = verts(nodes(triVerts(i,j)),1);
					triNodeIKY = verts(nodes(triVerts(i,j)),2);
					t = 1;
					while(t<=M)
						if((verts(nodes(t),1)==triNodeIKX)&&(verts(nodes(t),2)==triNodeIKY))
							break;
						endif
						t++;
					endwhile
					%Step 18
					if(l<=n)
						if(t<=n)
							alpha(l,t) = alpha(l,t) + J(i,k,j);
							alpha(t,l) = alpha(t,l) + J(i,k,j);
						else
							beta(l) = beta(l) - gamma(t)*J(i,k,j);
						endif
					else
						if(t<=n)
							beta(t) = beta(t)-gamma(l)*J(i,k,j);
						endif
					endif
					j++;
				endwhile
			endif
			%Step 19
			if(l<=n)
				alpha(l,l) = alpha(l,l) + J(i,k,k);
				beta(l) = beta(l) + I(i,k);
			endif
			k++;
		endwhile
		i++;
	endwhile
	disp(alpha);
	disp(beta');
	solution = alpha \ beta';
	disp(solution);
	i = 1;
	while(i<=M)
		j = 1;
		while(j<=3)
			printf("Triangle %i, a(%i): %4g, b(%i): %4g, c(%i): %4g\n",i,j,a(i,j),j,b(i,j),j,c(i,j));
			j++;
		endwhile
		i++;
	endwhile
	return;
endfunction