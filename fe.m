% There are ONLY Dirichlet boundary conditions in this algorithm,
% so boundary value problems that this algorithm can be used to solve are limited.
% The assembly of the final approximate model was done by hand using Mathematica
% If anyone knows of a way to do this deterministically in matlab, let me know!

%g - Dirichlet boundary condition
%verts - Vertices of the triangulation
%nodes - each node corresponds to the row of verts which contains the vertex for that node
%triVerts - each row has 3 values - the three node numbers of the triangle
%K - Interior triangles
%M - All triangles (total triangles - M)
%n - Nodes on the interior
%m - Total nodes
function f = fe(f,g,verts,nodes,triVerts,K,M,n,m)
	gamma = zeros(1,m);
	l = n+1;
	while(l<=m)
		gamma(l) = g(verts(nodes(l),1),verts(nodes(l),2));
		l++;
	endwhile
	alpha = zeros(n,n);
	beta = zeros(1,n);
	delta = zeros(1,M);
	ident = @(x,y)(1);
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
				firstTerm = b(i,j)*b(i,k)*area;
				secondTerm = c(i,j)*c(i,k)*area;
				z(i,j,k) = firstTerm+secondTerm;
				hInt = @(x,y)(f(x,y)*nj(x,y));
				H(i,j) = area*(hInt(midx1,midy1)+hInt(midx2,midy2)+hInt(midx3,midy3));
				k++;
			endwhile
			j++;
		endwhile
		i++;
	endwhile
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
	solution = alpha \ beta';
	disp(solution);
	disp(gamma);
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

%The following is a test of the finite element method to solve the steady-state
%heat equation with boundary conditions 0C on 3 sides and 100C on the fourth side.
%There are 16 triangles in the discretization.
f = @(x,y)(0);
function res = gp(x,y)
	if(x==0&&y<10)
		res = 0;
	elseif(x==10&&y<10)
		res = 0;
	elseif(y==0)
		res = 0;
	else
		res = 100;
	endif
endfunction
g = @(x,y)(gp(x,y));
K = 8;
M = 16;
m = 13;
n = 5;
verts = [2.5,2.5;5,5;7.5,2.5;7.5,7.5;2.5,7.5;0,10;0,5;0,0;5,0;10,0;10,5;10,10;5,10];
nodes = [1;2;3;4;5;6;7;8;9;10;11;12;13];
triVerts = [1,9,2;9,3,2;2,3,11;2,11,4;2,4,13;2,13,5;2,5,7;2,7,1;9,10,3;3,10,11;11,12,4;4,12,13;5,13,6;5,6,7;7,8,1;8,9,1];