function z = SplineCuadLocal(x, y, d_k, k)
	d = d_k;
	h = diff(x);
	p = diff(y)./h;

  #Recorremos todos los nodos de n+1 en adelante:
 if (k != (length(x) - 1))
  for i = (k+1):(length(x)-1)
		q = (p(i) - d)/h(i);
		v = [x(i) x(i)];
		s(i,:) = [0, 0, y(i)] + [0, d, -d*x(i)] +q*poly(v);
		d = polyval(polyder(s(i,:)),x(i+1));
  end
  d = d_k;
 endif

  #Recorremos todos los nodos desde n hasta el 1:
 if (k != 0)
  for j = k:-1:1
		q = (d - p(j))/h(j);
		v = [x(j) x(j+1)];
		s(j,:) = [0 0 y(j)] + [0 p(j) -p(j)*x(j)] + q*poly(v);
		d = polyval(polyder(s(j,:)), x(j));
  end
 endif

  for i = 1:(length(x) - 1)
    s(i,:) = polyaffine(s(i,:), [-x(i), 1]);
  end

  z = mkpp(x, s);
end
