function s = SplineNat(x, y)
  n = length(x) - 1;

  h = diff(x);     #h_i
  p = diff(y)./h;  #p_i

  twoes  = 2*ones(1,n-1);
  mu     = (h(1:end-1)./(h(1:end-1) + h(2:end)));
  lambda = ones(1,n-1) - mu;
  g      = diff(p)./(h(1:end-1) + h(2:end));

  if(n < 3)
    A = 2;
  else
    A = diag(twoes) + diag(lambda(1:end-1),1) + diag(mu(2:end),-1);
  end

  m = A \ g';
  s = ppint(ppint(SplineLineal(x, [0 m' 0])));
end
