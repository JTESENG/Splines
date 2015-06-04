function z = SplineCuadLocal(x, y, d_k, k)
  s = zeros(length(x)-1, 3);
  d = d_k;

  #Recorremos todos los nodos de n+1 en adelante:

  for i = (k+1):(length(x)-1)
		p = (y(i+1)-y(i))/(x(i+1)-x(i));
		q = (p-d)/(x(i+1)-x(i));
		v = [x(i) x(i)];
		s(i,:) = [0, 0, y(i)]+[0, d, -d*x(i)]+q*poly(v);
		d = polyval(polyder(s(i,:)),x(i+1));
	end
    d = d_k;

  #Recorremos todos los nodos desde n hasta el 1:

  for j = k:-1:1
		p = (y(j+1)-y(j))/(x(j+1)-x(j))
		q = (d-p)/(x(j+1)-x(j))
		v = [x(j) x(j+1)]
		s(j,:) = [0 0 y(j)]+[0 p -p*x(j)]+q*poly(v);
		d = polyval(polyder(s(j,:)), x(j))
  end

  for i = 1:length(s)
    s(i,:) = polyaffine(s(i,:), [-x(i), 1]);
  end

  z = mkpp(x, s);
end
