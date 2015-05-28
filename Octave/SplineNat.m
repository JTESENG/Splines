function s = SplineNat(x, y)
  n = length(x) - 1

  h = diff(x);  #h_i
  p = h./diff(y); #p_i

  u = h./(h + h(2:end)); # mu_i
  v = ones(1,n) - u;     # lambda_i
  g = diff(p)./(h + h(2:end))  # gamma_i

  A = diag(2*ones(1,n+1)) + diag(v,1) + diag(u,-1)
  m = A \ g';

  s = ppint(ppint(SplineLineal(x, m)));
end
