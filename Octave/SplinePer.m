function z = SplinePer (x, y) 
  n     = length(x);
  twoes = 2*ones(1,n);
  h     = diff(x);

  for i=1:(n-2) 
  	lambda(i+1) = h(i+1)/(h(i+1)+h(i));
  endfor

  for i=1:(n-2)
  	mu(i) = mu(i)+h(i)/(h(i+1)+h(i));
  endfor

  A = diag(mu,-1) + diag(twoes) + diag(lambda,1); 
  A(1, n) = h(n-1)/(h(n-1)+h(n-2));
  
  dd1 = diff(y)./diff(x); #Primeras derivadas

  #Segundas derivadas
  dd2 = zeros(n-2,1);
  for i=1:(n-2)
  	dd2(i) = (dd1(i+1)-dd1(i))/(x(i+2)-x(i));
  endfor
  
  A(1,:)=h(1)*A(1,:);   #Ajustes de la matriz
  A(1,1)=A(1,1)-h(1)/3;
  A(1,2)=A(1,2)+h(1)/6;
  A(n,:)=h(n-1)*A(1,:);
  A(n,1)=A(1,1)+h(1)/3;
  A(n,2)=A(1,2)-h(1)/6;

  gamma(1)=(h(1)-dd1(1))/h(1) + (dd1(1)-h(n))/h(n);
  for i=1:(n-2)
  	gamma(i+1)=6*dd2(i);
  endfor

  m=A\gamma';

	for i = 1:n
	  p  = (-m(i)/6)  * poly([x(i+1), x(i+1), x(i+1)]);
	  p += (m(i+1)/6) * poly([x(i), x(i), x(i)]);
	  p += (y(i)-(m(i)*h(i).^2/6))*[0 0 -1 x(i+1)];
	  p += (y(i+1)-(m(i+1)*h(i).^2/6))*[0 0 1 -x(i)];
	  p *= 1/h(i);
	  B(i,:) = polyaffine(p, [-x(i) 1]);
  	endfor
  	s = mkpp(x, B);
endfunction
