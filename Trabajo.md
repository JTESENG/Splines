---
title: Ampliación de interpolación con Splines
author: [Miguel Anguita Ruiz, Pablo Baeyens Fernández, Pablo David Medina Sánchez, Ruben Morales Pérez, Francisco Javier Morales Piqueras]
lang: spanish
header-includes:
toc: false
numbersections: false
fontsize: 11pt
geometry: margin=1in
---

# Splines cuadráticos

El espacio de splines de clase 2 se nota $S_2(x_1,x_2,...,x_n)$ Una base sería ${1,x,x^2,(x-x_1)_+,(x-x_2)_+,...,(x-x_n)_+}$

## Descripción del espacio de splines cuadráticos
El espacio de splines de clase 2 con n nodos se denota $S_2(x_1,x_2,...,x_n)$. Los splines de clase 2 están constituídos por parábolas de forma que además de tener una función continua, su derivada también lo es, por lo tanto, para $i=1,...,n-1$ tenemos la siguiente condición:

\begin{tabular}{|c|}
   \hline
   $s_i(x_i)=s_{i+1}(x_{i+1})$\\
   $s_i'(x_i)=s_{i+1}'(x_{i+1})$\\
   \hline
\end{tabular}

Proposición:<!---Ponerlo esto mejor, con mejor formato-->
El conjunto $S_2(x_1,x_2,...,x_n)$ satisface las propiedades siguienes:
\begin{enumerate}
  \item Es un espacio vectorial con $\dim(S_2(x_1,x_2,...,x_n))$
\end{enumerate}


## Ejemplos


# Splines cúbicos

## Construcción a partir de los valores de $s''(x)$ en los nodos $\{x_i\}$

## Propiedades de minimización

## Ejemplos

# Implementación en ordenador: Octave

## Splines cúbicos

## Splines cuadráticos

### Resolviendo a trozos

### Mediante sistema de ecuaciones

Utilizando el sistema que vimos anteriormente, podemos definir fácilmente una
función que calcule los coeficientes de un spline cuadrático de clase 1:

```octave
function s = coefsSpline(x, y, d_k, k)
  # Número de intervalos
  n = length(x) - 1;

  # 1, x, x²
  A(:,1) = [ones(n+1,1); 0];
  A(:,2) = [x'         ; 1];
  A(:,3) = [x'.^2      ; 2.*x(k+1)];

  # Potencias truncadas
  for j = 4 : n + 2
    pot    =  @(t) (t > x(j-2)) .* (t - x(j-2));
    A(:,j) = [pot(x').^2; 2.*pot(x(k+1))];
  end

  # Resolución del sistema
  s = A \ [y' ; d_k];

end
```

Definimos la matriz columna a columna: las **3 primeras columnas** corresponden a
los valores de `x` en $1, x, x^2$, así como los valores de la derivada:

```octave
A(:,1) = [ones(n+1,1); 0];
A(:,2) = [x'         ; 1];
A(:,3) = [x'.^2      ; 2.*x(k+1)];
```

Una vez hecho esto pasamos a las **potencias truncadas**. Para ello, en cada
columna:

- Definimos una función `pot` correspondiente a la potencia truncada en el valor
correspondiente: `pot    =  @(t) (t > x(j-2)) .* (t - x(j-2))`. Como *Octave*
tiene tipos dinámicos convertirá `(t > x(j-2))` a $1$ o $0$. De esta forma,
$pot(x) = (x - x_{j-1})_{+}$.

- Aplicamos `pot` a `x` en cada columna, añadiendo el valor de la derivada.

Definida la matriz del sistema podemos calcular finalmente la solución empleando
la **división izquierda**.
