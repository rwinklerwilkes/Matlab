%Linear Shooting Using 4th Order RK
function w = ls(p,q,r,a,b,alpha,beta,N)
	u = zeros(2,N+1);
	v = zeros(2,N+1);
	w = zeros(2,N);
	k = zeros(4,2);
	kk = zeros(4,2);
	h = (b-a)/N;
	u(1,1) = alpha;
	u(2,1) = 0;
	v(1,1) = 0;
	v(2,1) = 1;
	i = 1;
	while(i<=N)
		x = a+(i-1)*h;
		k(1,1) = h*u(2,i);
		k(1,2) = h*(p(x)*u(2,i)+q(x)*u(1,i)+r(x));
		k(2,1) = h*(u(2,i)+.5*k(1,2));
		k(2,2) = h*(p(x+h/2)*(u(2,i)+.5*k(1,2)) + q(x+h/2)*(u(1,i)+.5*k(1,1))+r(x+h/2));
		k(3,1) = h*(u(2,i) + .5*k(2,2));
		k(3,2) = h*(p(x+h/2)*(u(2,i)+.5*k(2,2)) + q(x+h/2)*(u(1,i) + .5*k(2,1)) + r(x+h/2));
		k(4,1) = h*(u(2,i) + k(3,2));
		k(4,2) = h*(p(x+h)*(u(2,i) + k(3,2))+q(x+h)*(u(1,i) + k(3,1)) + r(x+h));
		u(1,i+1) = u(1,i) + (1/6)*(k(1,1) + 2*k(2,1) +2*k(3,1) +k(4,1));
		u(2,i+1) = u(2,i) + (1/6)*(k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2));
		kk(1,1) = h*v(2,i);
		kk(1,2) = h*(p(x)*v(2,i)+q(x)*v(1,i));
		kk(2,1) = h*(v(2,i)+.5*kk(1,2));
		kk(2,2) = h*(p(x+h/2)*(v(2,i)+.5*kk(1,2)) + q(x+h/2)*(v(1,i)+.5*kk(1,1)));
		kk(3,1) = h*(v(2,i) + .5*kk(2,2));
		kk(3,2) = h*(p(x+h/2)*(v(2,i)+.5*kk(2,2)) + q(x+h/2)*(v(1,i) + .5*kk(2,1)));
		kk(4,1) = h*(v(2,i) + kk(3,2));
		kk(4,2) = h*(p(x+h)*(v(2,i) + kk(3,2))+q(x+h)*(v(1,i) + kk(3,1)));
		v(1,i+1) = v(1,i) + (1/6)*(kk(1,1) + 2*kk(2,1) +2*kk(3,1) +kk(4,1));
		v(2,i+1) = v(2,i) + (1/6)*(kk(1,2) + 2*kk(2,2) + 2*kk(3,2) + kk(4,2));
		i++;
	endwhile
	w(1,1) = alpha;
	w(2,1) = (beta-u(1,N+1))/v(1,N+1);
	i = 1;
	while(i<=N+1)
		W1 = u(1,i) + w(2,1)*v(1,i);
		W2 = u(2,i) + w(2,1)*v(2,i);
		x = a+(i-1)*h;
		i++;
		printf("x: %4g, W1: %6g, W2: %6g\n",x,W1,W2);
	endwhile
	return;
endfunction