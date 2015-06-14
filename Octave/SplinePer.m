function s = SplinePer (x, y)
  n     = length(x);
  twoes = 2*ones(1,n);
  h     = diff(x);

  for i=1:(n-2)
  	lambda(i+1) = h(i+1)/(h(i+1)+h(i));
  endfor

  for i=1:(n-2)
  	mu(i) = h(i)/(h(i+1)+h(i));
  endfor

  A = diag([mu 0],-1) + diag(twoes) + diag(lambda,1);
  A(1,:) = [1 zeros(1,n-2) -1];
  A(n,:) = [-2*h(1) -h(1) zeros(1,n-4) -h(n-1) -2*h(n-1)];

  dd1 = diff(y)./h; #Primeras derivadas

  for i=1:(n-2)
  	dd2(i) = (dd1(i+1)-dd1(i))/(x(i+2)-x(i));
  endfor

  gamma = [0 6*dd2 -6*(dd1(1) - dd1(n-1))];
  m     = A\gamma';

	for i = 1:n-1
	  p  = (-m(i)/6)  * poly([x(i+1), x(i+1), x(i+1)]);
	  p += (m(i+1)/6) * poly([x(i), x(i), x(i)]);
	  p += (y(i)-(m(i)*h(i).^2/6))*[0 0 -1 x(i+1)];
	  p += (y(i+1)-(m(i+1)*h(i).^2/6))*[0 0 1 -x(i)];
	  p *= 1/h(i);
	  B(i,:) = polyaffine(p, [-x(i) 1]);
  	endfor
  	s = mkpp(x, B);
endfunction
