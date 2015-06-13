function s = SplinePer (x, y)
  n     = length(x);
  twoes = 2*ones(1,n);
  h     = diff(x);

  lambda= zeros(1,n-1);
  for i=1:(n-2)
  	lambda(i+1) = h(i+1)/(h(i+1)+h(i));
  endfor
  mu= zeros(1,n-1);
  for i=1:(n-2)
  	mu(i) = mu(i)+h(i)/(h(i+1)+h(i));
  endfor


  prim_fila=zeros(1,n);
  prim_fila(1)=1;
  prim_fila(n)=-1;
  ult_fila=zeros(1,n);
  ult_fila(1)=3*h(1);
  ult_fila(2)=-h(1);
  ult_fila(n-1)=-h(n-1);
  ult_fila(n)=-5*h(n-1);

  A = diag(mu,-1) + diag(twoes) + diag(lambda,1);

  A(1,:)=prim_fila;
  A(n,:)=ult_fila;
  dd1 = diff(y)./diff(x); #Primeras derivadas
  dd2 = zeros(n-2,1);     #Segundas derivadas
  for i=1:(n-2)
  	dd2(i) = (dd1(i+1)-dd1(i))/(x(i+2)-x(i));
  endfor

  gamma = zeros(1,n);
  for i=1:(n-2)
  	gamma(i+1)=6*dd2(i);
  endfor
  gamma(n)= -6*( (y(2)-y(1))/(x(2)-x(1)) - (y(n)-y(n-1))/(x(n)-x(n-1)) );

  m=A\gamma';

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
