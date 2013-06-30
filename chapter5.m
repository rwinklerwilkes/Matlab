#Euler's Method
function e = euler(f,a,b,N,init)
	h = (b-a)/N;
	t = a;
	e(1) = init;
	i = 1;
	while(i<=N)
		f(t,e(i));
		e(i+1) = e(i) + h*f(t,e(i));
		t = a+i*h;
		i++;
	endwhile
	return;
endfunction

#Midpoint Method (O(h^2))
function w = midpoint(f,t,tend,alpha,h)
	N = (tend-t)/h;
	w(1) = alpha;
	i = 1;
	while(i<=N)
		a = t+(i-1)*h;
		w(i+1) = w(i) + h*f((a + h/2),(w(i) + (h/2)*f(a,w(i))));
		i++;
	endwhile
	return;
endfunction

#Modified Euler Method (O(h^2))
function w = meuler(f,t,tend,alpha,h)
	N = (tend-t)/h;
	w(1) = alpha;
	i = 1;
	a = t+(i-1)*h;
	while(i<=N)
		aa = t+i*h;
		w(i+1) = w(i)+(h/2)*(f(a,w(i))+f(aa,w(i) + h*f(a,w(i))));
		a = aa;
		i++;
	endwhile
	return;
endfunction

#Heun's Method (O(h^3))
function w = heun(f,t,tend,alpha,h)
	N = (tend-t)/h;
	w(1) = alpha;
	i = 1;
	while(i<=N)
		a = t+(i-1)*h;
		w(i+1) = w(i)+(h/4)*(f(a,w(i))+3*f(a+(2*h)/3,w(i) + (2*h)/3*f(a+h/3,w(i)+(h/3)*f(a,w(i)))));
		i++;
	endwhile
	return;
endfunction

#Runge-Kutta Method (O(h^4))
function w = rk(f,t,tend,alpha,h)
	N = int32((tend-t)/h);
	w(1) = alpha;
	i = 1;
	a = t+(i-1)*h;
	while(i<=N)
		aNext= t+i*h;
		b = w(i);
		kOne = h*f(a,b);
		kTwo = h*f(a+h/2, b+kOne/2);
		kThree = h*f(a+h/2, b+kTwo/2);
		kFour = h*f(aNext,w(i)+kThree);
		w(i+1) = w(i) + (1/6)*(kOne + 2*kTwo + 2*kThree + kFour);
		a = aNext;
		i++;
	endwhile
	return;
endfunction

#Runge-Kutta Fehlberg Method (O(h^4))
function w = rkf(f,a,b,alpha,tol,hmin,hmax)
	t = a;
	w = alpha;
	h = hmax;
	FLAG = 1;
	disp(t),disp(w);
	while(FLAG == 1)
		k1 = h*f(t,w);
		k2 = h*f(t+(1/4)*h,w+(1/4)*k1);
		k3 = h*f(t+(3/8)*h,w+(3/32)*k1+(9/32)*k2);
		k4 = h*f(t+(12/13)*h,w+(1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
		k5 = h*f(t+h,w+(439/216)*k1 -8*k2 + (3680/513)*k3 - (845/4104)*k4);
		k6 = h*f(t+h/2,w-(8/27)*k1+2*k2-(3544/2565)*k3+(1859/4104)*k4-(11/40)*k5);
		R = (1/h)*abs((k1/360-(128/4275)*k3-(2197/75240)*k4+k5/50+(2*k6)/55));
		if(R<=tol)
			t = t+h;
			w = w + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - k5/5;
			disp("t, w, h"),disp(t),disp(w),disp(h);
		endif
		delta = .84*(tol/R)**.25;
		if(delta<=.1)
			h = .1*h;
		elseif(delta>=4)
			h = 4*h;
		else
			h = delta*h;
		endif
		if(h>hmax)
			h = hmax;
		endif;
		if(t>=b)
			FLAG = 0;
		elseif(t+h > b)
			h = b-t;
		elseif(h < hmin)
			FLAG = 0;
			disp("Minimum h value exceeded");
		endif
	endwhile
endfunction

#Predictor-Corrector Method
function p = pc(f,a,b,N,alpha)
	h = (b-a)/N;
	t(1) = a;
	t(2) = a+h;
	t(3) = a+2*h;
	t(4) = a+3*h;
	p = rk(f,a,a+3*h,alpha,h);
	i = 5;
	while(i<=(N+1))
		t(i) = a + (i-1)*h;
		#Predict next value
		p(i) = p(i-1) + (h/24)*(55*f(t(i-1),p(i-1)) - 59*f(t(i-2),p(i-2)) + 37*f(t(i-3),p(i-3)) - 9*f(t(i-4),p(i-4)));
		#Correct using prediction
		p(i) = p(i-1) + (h/24)*(9*f(t(i),p(i)) + 19*f(t(i-1),p(i-1)) -5*f(t(i-2),p(i-2)) +f(t(i-3),p(i-3)));
		i++;
	endwhile
	return;
endfunction

#Modified Runge-Kutta for use in Adaptive Predictor-Corrector
function r = prk(f,h,xVals,yVals)
	j = 2;
	while(j<=4)
		kOne = h*f(xVals(j-1),yVals(j-1));
		kTwo = h*f(xVals(j-1)+h/2,yVals(j-1)+kOne/2);
		kThree = h*f(xVals(j-1)+h/2,yVals(j-1)+kTwo/2);
		kFour = h*f(xVals(j-1)+h,yVals(j-1)+kThree);
		yVals(j) = yVals(j-1)+(kOne+2*kTwo+2*kThree+kFour)/6;
		xVals(j) = xVals(1) + (j-1)*h;
		j++;
	endwhile
	r = vertcat(xVals,yVals);
	return;
endfunction

#This needs work, not yet finished
#Adaptive Predictor-Corrector Method
function p = apc(f,a,b,alpha,tolerance,minStep,maxStep)
	h = maxStep;
	t(1) = a;
	t(2) = a+h;
	t(3) = a+2*h;
	t(4) = a+3*h;
	p = prk(f,h,t,alpha);
	#PRK will be a 2xN matrix with first row t values and second row W values
	i = 5;
	FLAG = 1;
	LAST = 0;
	NFLAG = 1;
	while(FLAG==1)
		t = p(1,i-1) + h;
		WP = p(2,i-1) + (h/24)*(55*f(p(1,i-1),p(2,i-1)) - 59*f(p(1,i-2),p(2,i-2)) + 37*f(p(1,i-3),p(2,i-3)) - 9*f(p(1,i-4),p(2,i-4)));
		WC = p(2,i-1) + (h/24)*(9*f(t,WP) + 19*f(p(1,i-1),p(2,i-1)) -5*f(p(1,i-2),p(2,i-2)) +f(p(1,i-3),p(2,i-3)));
		delta = 19*abs(WC-WP)/(270*h);
		if(delta<=tolerance)
			p(1,i) = t;
			p(2,i) = WC;
			if(NFLAG==1)
				j = i-3;
				while(j<=i)
					printf("j: %6d t: %6.6g w: %6.6g h: %6.6g\n",j,p(1,j),p(2,j),h);
					j++;
				endwhile
			else
				printf("j: %6d t: %6.6g w: %6.6g h: %6.6g\n",j,p(1,j),p(2,j),h);
			endif
			if(LAST == 1)
				FLAG=0;
			else
				i++;
				NFLAG = 0;
				if((delta < .01*tolerance)|| p(1,i-1)+h > b)
					q = (tolerance/(2*delta))**.25;
					if(q>4)
						h = 4*h;
					else
						h = q*h;
					endif
					if(h > maxStep)
						h = maxStep;
					endif
					if(p(1,i-1)+h > b)
						h = (b-p(1,i-1))/4;
						LAST = 1;
					endif
					p(1:2,(i-1):(i+2)) = [prk, prk(f,h,p(1,(i-1):(i+2)),p(2,(i-1):(i+2)))(1:2,2:4)];
					NFLAG = 1;
					i = i+3;
				endif
			endif
		else
			q = (tolerance/(2*delta))**.25
			if(q<.1)
				h = .1*h;
			else
				h = q*h;
			endif
			if(h<minStep)
				FLAG = 0;
				disp("Minimum step size exceeded");
			else
				if(NFLAG == 1)
					i = i-3;
					p(1:2,(i-1):(i+2)) = prk(f,h,p(1,(i-1):(i+2)),p(2,(i-1):(i+2)));
					i = i+3;
					NFLAG = 1;
				endif
			endif
		endif
	endwhile
	return;
endfunction

#Extrapolation Method
function e = extrapolate(f,a,b,alpha,tolerance,hmin,hmax)
	NK = [2,4,6,8,12,16,24,32];
	TO = a;
	WO = alpha;
	h = hmax;
	FLAG = 1;
	i = 1;
	while(i<=7)
		j = 1;
		while(j<=i)
			Q(i,j) = (NK(i+1)/NK(j))**2;
			j++;
		endwhile
		i++;
	endwhile
	while(FLAG==1)
		k = 1;
		NFLAG = 0;
		while(k<=4&&NFLAG==0)
			HK = h/NK(k);
			T = TO;
			W2 = WO;
			W3 = W2+HK*f(T,W2);
			T = TO + HK;
			j = 1;
			while(j<=(NK(k)-1))
				W1 = W2;
				W2 = W3;
				W3 = W1+2*HK*f(T,W2);
				T = TO + (j+1)*HK;
				j++;
			endwhile
			y(k) = (W3 + W2 + HK*f(T,W3))/2;
			if(k>=2)
				j = k;
				v = y(1);
				while(j>=2)
					y(j-1) = y(j) + (y(j)-y(j-1))/(Q(k-1,j-1)-1);
					j--;
				endwhile
				if(abs(y(1) - v)<=tolerance)
					NFLAG = 1;
				endif
			endif
			k++;
		endwhile
		k--;
		if(NFLAG==0)
			h = h/2;
			if(h<hmin)
				disp("Minimum step size exceeded");
				FLAG = 0;
			endif
		else
			WO = y(1);
			TO = TO+h;
			printf("t: %6.6g, w: %12.10g, h:%10.10g\n",TO,WO,h);
			if(TO >= b)
				FLAG = 0;
			elseif(TO+h > b)
				h = b-TO;
			elseif(k<=3 && h<.5*hmax)
				h = 2*h;
			endif
		endif
	endwhile
endfunction

#Runge-Kutta Method for Systems of Differential Equations
function r = rks(a,b,fOne,fTwo,alpha,N)
	r = alpha;
	h = (b-a)/N;
	t = a;
	for j = 1:2
		w(j) = alpha(j);
	endfor
	for i = 1:N
		for j = 1:2
			if(j==1)
				k(1,j) = h*fOne(t,w(1),w(2));
			else
				k(1,j) = h*fTwo(t,w(1),w(2));
			endif
		endfor
		for j = 1:2
			if(j==1)
				k(2,j) = h*fOne(t+h/2,w(1)+.5*k(1,1),w(2)+.5*k(1,2));
			else
				k(2,j) = h*fTwo(t+h/2,w(1)+.5*k(1,1),w(2)+.5*k(1,2));
			endif
		endfor
		for j = 1:2
			if(j==1)
				k(3,j) = h*fOne(t+h/2,w(1)+.5*k(2,1),w(2)+.5*k(2,2));
			else
				k(3,j) = h*fTwo(t+h/2,w(1)+.5*k(2,1),w(2)+.5*k(2,2));
			endif
		endfor
		for j = 1:2
			if(j==1)
				k(4,j) = h*fOne(t+h,w(1)+k(3,1),w(2)+k(3,2));
			else
				k(4,j) = h*fTwo(t+h,w(1)+k(3,1),w(2)+k(3,2));
			endif
		endfor
		for j = 1:2
			w(j) = w(j) + (k(1,j) + 2*k(2,j) + 2*k(3,j) + k(4,j))/6;
		endfor
		t = a+i*h;
		printf("t: %10.10g w1: %12.10g w2: %12.10g\n",t,w(1),w(2));
		r = vertcat(r,[w(1),w(2)]);
	endfor
	return;
endfunction

function r = rks3(a,b,fOne,fTwo,fThree,alpha,N)
h = (b-a)/N;
	t = a;
	for j = 1:3
		w(j) = alpha(j);
	endfor
	for i = 1:N
		for j = 1:3
			if(j==1)
				k(1,j) = h*fOne(t,w(1),w(2),w(3));
			elseif(j==2)
				k(1,j) = h*fTwo(t,w(1),w(2),w(3));
			else
				k(1,j) = h*fThree(t,w(1),w(2),w(3));
			endif
		endfor
		for j = 1:3
			if(j==1)
				k(2,j) = h*fOne(t+h/2,w(1)+.5*k(1,1),w(2)+.5*k(1,2),w(3)+.5*k(1,3));
			elseif(j==2)
				k(2,j) = h*fTwo(t+h/2,w(1)+.5*k(1,1),w(2)+.5*k(1,2),w(3)+.5*k(1,3));
			else
				k(2,j) = h*fThree(t+h/2,w(1)+.5*k(1,1),w(2)+.5*k(1,2),w(3)+.5*k(1,3));
			endif
		endfor
		for j = 1:3
			if(j==1)
				k(3,j) = h*fOne(t+h/2,w(1)+.5*k(2,1),w(2)+.5*k(2,2),w(3)+.5*k(2,3));
			elseif(j==2)
				k(3,j) = h*fTwo(t+h/2,w(1)+.5*k(2,1),w(2)+.5*k(2,2),w(3)+.5*k(2,3));
			else
				k(3,j) = h*fThree(t+h/2,w(1)+.5*k(2,1),w(2)+.5*k(2,2),w(3)+.5*k(2,3));
			endif
		endfor
		for j = 1:3
			if(j==1)
				k(4,j) = h*fOne(t+h,w(1)+k(3,1),w(2)+k(3,2),w(3)+k(3,3));
			elseif(j==2)
				k(4,j) = h*fTwo(t+h,w(1)+k(3,1),w(2)+k(3,2),w(3)+k(3,3));
			else
				k(4,j) = h*fThree(t+h,w(1)+k(3,1),w(2)+k(3,2),w(3)+k(3,3));
			endif
		endfor
		for j = 1:3
			w(j) = w(j) + (k(1,j) + 2*k(2,j) + 2*k(3,j) + k(4,j))/6;
		endfor
		t = a+i*h;
		printf("t: %10.10g w1: %12.10g w2: %12.10g w3: %12.10g\n",t,w(1),w(2),w(3));
	endfor
endfunction

#Implicit Trapezoidal Rule with Newton Iteration
#Used for solving stiff differential equations
function tr = tni(f,fy,a,b,N,alpha,tol,maxIter)
	h = (b-a)/N;
	w(1) = alpha;
	t(1) = a;
	i = 1;
	while(i<=N)
		k1 = w(i) + (h/2)*f(t(i),w(i));
		ww = k1;
		j = 1;
		FLAG = 0;
		while(FLAG==0)
			w(i+1) = ww - (ww - (h/2)*f(t(i)+h,ww)-k1)/(1-(h/2)*fy(t(i)+h,ww));
			if(abs(w(i+1)-ww) < tol)
				FLAG = 1;
			else
				j++;
				ww = w(i+1);
				if(j>maxIter)
					disp("Maximum number of iterations exceeded");
					return;
				endif
			endif
		endwhile
		t(i+1) = a+i*h;
		i++;
	endwhile
	tr = vertcat(t,w);
	return;
endfunction