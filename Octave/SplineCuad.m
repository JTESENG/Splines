function s = SplineCuad(x, y, d_k, k)
  # Número de intervalos
  n = length(x) - 1;

  # 1, x, x²
  A(:,1) = [ones(n+1,1); 0];
  A(:,2) = [x'         ; 1];
  A(:,3) = [x'.^2      ; 2.*x(k+1)];

  # Potencias truncadas
  for j = 4 : n + 2
    t       =  @(s) (s > x(j-2)) .* (s - x(j-2));
    A(:, j) = [t(x').^2; 2.*t(x(k+1))];
  end

  # Resolución del sistema
  sol = A \ [y' ; d_k];

  for k = 1:n
    p = sol(3:-1:1)';

    for l = 2:k
      p += sol(l+2).*[1, -2.*x(l), x(l).^2];
    end

    B(k, :) = polyaffine(p,[-x(k) 1]);
  end

  s = mkpp(x,B);
end
