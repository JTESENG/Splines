function s = SplineNat(x, y)
  n = length(x) - 1;

  h = diff(x);     #h_i
  p = diff(y)./h;  #p_i

  twoes  = 2*ones(1,n-1);
  mu     = (h(1:end-1)./(h(1:end-1) + h(2:end)));
  lambda = ones(1,n-1) - mu;
  gamma  = diff(p)./(h(1:end-1) + h(2:end));

  if(n < 3)
    A = 2;
  else
    A = diag(twoes) + diag(lambda(1:end-1),1) + diag(mu(2:end),-1)
  end

  m = [0, (A \ (6.*gamma'))', 0];

  for i = 1:n
    p  = (-m(i)/6)  * poly([x(i+1), x(i+1), x(i+1)]);
    p += (m(i+1)/6) * poly([x(i), x(i), x(i)]);
    p += (y(i)-(m(i)*h(i).^2/6))*[0 0 -1 x(i+1)];
    p += (y(i+1)-(m(i+1)*h(i).^2/6))*[0 0 1 -x(i)];
    p *= 1/h(i);
    B(i,:) = polyaffine(p, [-x(i) 1]);
  end
  s = mkpp(x, B);
end
