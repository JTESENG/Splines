---
title: Ampliación de interpolación con Splines
author: [Miguel Anguita Ruiz, Pablo Baeyens Fernández, Pablo David Medina Sánchez, Ruben Morales Pérez, Francisco Javier Morales Piqueras]
lang: spanish
header-includes:
	\usepackage{mathrsfs}
	\usepackage{amsthm}
	\newtheorem*{proposition}{Proposición}
	\newtheorem*{teorema}{Teorema}
	\theoremstyle{definition}
	\newtheorem*{definicion}{Definición}
	\newtheorem*{problema}{Problema}
toc: true
numbersections: true
fontsize: 11pt
geometry: margin=1in
---

\pagebreak

# Splines cuadráticos

El espacio de splines de clase 2 se nota $S_2(x_1,x_2,...,x_n)$.
Una base sería $\{1,x,x^2,(x-x_1)_+,(x-x_2)_+,...,(x-x_n)_+\}$

## Descripción del espacio de splines cuadráticos
El espacio de splines de clase 2 con $n$ nodos se denota $S_2(x_1,x_2,...,x_n)$.
Los splines de clase 2 están constituidos por parábolas, de forma que además de
tener una función continua, su derivada también lo es.
Por lo tanto, para $i=1,...,n-1$ tenemos la siguiente condición:

\begin{tabular}{|c|}
   \hline
   $s_i(x_i)=s_{i+1}(x_{i+1})$\\
   $s_i'(x_i)=s_{i+1}'(x_{i+1})$\\
   \hline
\end{tabular}

\begin{proposition}
El conjunto $S_2(x_1,x_2,...,x_n)$ satisface las propiedades siguienes:
\begin{enumerate}
  \item Es un espacio vectorial con $\dim(S_2(x_1,x_2,...,x_n))$
\end{enumerate}
\end{proposition}

## Interpolación con splines cuadráticos

## Error en los splines cuadráticos

\begin{teorema}
Sean $f \in C^2([a,b])$, $\{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$, $s \in S_2^1(\{x_i\}_{i = 0...n})$ spline para $f$, 

$h = max\{x_i - x_{i-1}\}_{i = 1...n}$, $E = f - s$. Además, sea $M >0$ tal que:

\[ M \geq Sup \{|f''(x) - f''(y)| \; : \; |x - y| \leq h, \; x,y \in [a,b] \}\]

Entonces, se verifica, para todo $x \in [a,b]$:

\begin{equation}
E(x) \leq \frac{h^2M}{2}
\end{equation}

\end{teorema}

La demostración, así como cotas para las derivadas y cotas más precisas en función
de la localización de $x$ puede encontrarse en *Quadratic Interpolatory Splines*,
W. Kammerer, G. Reddien y R.S. Varga, (1973).

## Ejemplos

\pagebreak

# Splines cúbicos

## Construcción a partir de los valores de $s''(x)$ en los nodos $\{x_i\}$

## Propiedades de minimización

Comenzamos planteando un problema de minimización sobre $(C^2([a,b]), || \cdot ||)$, con la norma definida de la forma usual:

\begin{equation}
|| f || = \sqrt{ \int_a^b f(x)^2 dx }
\end{equation}

El problema es aproximar una función de clase 2 con funciones que la interpolen en
unos nodos y cuyas derivadas en los extremos coincidan:
 
\begin{problema}
Sea $f \in C^2([a,b])$, $P \in \mathscr{P}([a,b])$. Sea $H \subset C^2([a,b])$ definido por:
\[H = \{g \in C^2([a,b]) \; : \; \forall p \in P \; g(p) = f(p) \text{ y } \; g'(a) = f'(a), \; g'(b) = f'(b)\} \]

Hallar $u \in H$ tal que $||f - u||$ sea mínimo.
\end{problema}


<!--\begin{teorema}[Minimización]
Sea $f \in C^2([a,b])$, $P \in \mathscr{P}([a,b])$. Se verifica:
\begin{equation}
\int_a^b f(x)^2 dx \geq \int_a^b f(x)^2 dx
\end{equation}
\end{teorema}-->



### Cota de error en los splines cúbicos

\begin{teorema}
Sea $f \in C^4([a,b])$, $n \in \mathbb{N}$, $P = \{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$ y  $s \in S_3^1(P)$ spline para $f$. Además, sean $h = max\{x_i - x_{i-1}\}_{i = 1...n}$, $M > 0$ cota superior de $|f^{iv)}|$ en $[a,b]$ y $E = f - s$, $x \in [a,b]$.

Se verifica:

\begin{equation} \label{eq:errorS31}
|E(x)| \leq \frac{5M}{384}h^4
\end{equation}

\end{teorema}

La demostración, así como cotas para las derivadas, puede consultarse en *Optimal Error Bounds for Cubic Spline Interpolation*, Charles Hall y Weston Meyer, (1976).

## Ejemplos

\pagebreak

# Implementación en ordenador: Octave

Hemos implementado las siguientes funciones en Octave:

1. `SplineLineal`: Calcula spline **lineal**. *(Usado en los splines cúbicos)*
1. `Spline31` : Calcula spline de **clase 1**.
1. `SplineNat`: Calcula spline **natural**.
2. `SplinePer`: Calcula spline **periódico**.
3. `SplineSuj`: Calcula spline **sujeto**.
8. `SplineCuad`: Calcula spline **cuadrático** de clase 1.

## Spline Lineal

La función que nos permite calcular un spline lineal es muy 
```octave
function s = SplineLineal(x,y)
  p = diff(y)./diff(x);
  A = [p' y(1:end-1)'];
  s = mkpp(x,A);
end
```

## Splines cuadráticos

Utilizando el sistema que vimos anteriormente, podemos definir fácilmente una
función que calcule los coeficientes de un spline cuadrático de clase 1:

```octave
function s = coefsSplineCuad(x, y, d_k, k)
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

\pagebreak
\appendix

# Definiciones y notación

\begin{definicion}
Sea $I \subset \mathbb{R}$ un intervalo cerrado y acotado con extremos $a,b$:

\begin{itemize}
\item Una \textbf{partición} $P$ de $I$ es un subconjunto finito de $I$ con $a,b\in P$ .
\item $\mathscr{P}(I)$ es el conjunto de todas las particiones de $I$.
\end{itemize}
\end{definicion}

\vspace*{3pt}

\begin{definicion} Sea $a \in \mathbb{R}$, $n \in \mathbb{N}$. La 
\textbf{potencia truncada} en $a$ de grado $n$, $(x - a)_+^n$ viene dada por: 

\[
 (x - a)_+^n =
  \begin{cases}
      0 			& \text{si } x \leq a \\
   (x - a)^n   & \text{si } x > a
  \end{cases}
\]
\end{definicion}

Cualquier potencia truncada de grado $n$ es de clase $n - 1$, y su derivada de 
orden $n$ presenta una discontinuidad en $a$. La derivada de $(x - a)_+^n$ en
$x$ es $n(x - a)_+^{n-1}$.

Su implementación en Octave es bastante sencilla: dados `a` y `n`, podemos definir
la potencia truncada como función anónima de la siguiente forma:

```octave
pot = @(x) (x > a) * (x - a)^n
```

Como Octave tiene tipos dinámicos convertirá `(x > a)` a $1$ si $x > a$ y a $0$
en otro caso.





