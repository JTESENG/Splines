function z = SplineCuadLocal(x, y, d_k, k)
	s = zeros(length(x)-1, 3);
  d = d_k;

  #Recorremos todos los nodos de n+1 en adelante:

  for i = (k+1):length(x)
		p = (y(i)-y(i-1))/(x(i)-x(i-1));
		q = (p-d)/(x(i)-x(i-1));
		v = [x(i-1) x(i-1)];
		s(i-1,:) = [0 0 y(i-1)]+[0 d -d*x(i-1)]+q*poly(v);
		d = 2*p-d;
	end
    d = d_k;

  #Recorremos todos los nodos desde n hasta el 1:

  for i = 0:(k-2)
		j = k-i;
		p = (y(j)-y(j-1))/(x(j)-x(j-1));
		q = (d-p)/(x(j)-x(j-1));
		v = [x(j-1) x(j)];
		s(j-1,:) = [0 0 y(j-1)]+[0 p -p*x(j-1)]+q*poly(v);
  end

  for i = 1:length(s)
    s(i,:) = polyaffine(s(i,:), [-x(i), 1]);
  end
    z = mkpp(x, s);
end
