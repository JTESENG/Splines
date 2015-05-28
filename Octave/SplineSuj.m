function s = SplineSuj (x, y, d_1, d_n)
  n     = length(x) - 1;
  twoes = 2*ones(1,n+1);
  h     = diff(x);

  mu     = [(h(1:end-1)./(h(1:end-1) + h(2:end))) 1];
  lambda = ones(1,n) - [0 mu(1:end-1)];

  A = diag(mu,-1) + diag(twoes,0) + diag(lambda,1);

  dd1 = diff(y)./h;
  for i=1:(n-1)
          dd2(i)=(dd1(i+1)-dd1(i))/(x(i+2)-x(i));
  end

  gamma(1)   = 6*(dd2(1)-d_1)/(x(2)-x(1));
  gamma(2:n) = 6*dd2(1:n-1);
  gamma(n+1) = 6*(d_n-dd2(n-1))/(x(n+1)-x(n));

  m = A \ gamma';
  s = ppint(ppint(SplineLineal(x, m')));
end
