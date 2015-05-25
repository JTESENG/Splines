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

##Introducción a los splines

Un spline es una curva diferenciable definida a trozos en un intervalo [a,b]. Suponiendo que tenemos este intervalo, se define una partición P del intervalo anterior como P = {X0=a < X1 <...< Xn-1 < Xn=b}. Esta es la base de todo spline y el punto de partida para definir nuestro área de estudio: el spline cuadrático.

## Descripción del espacio de splines cuadráticos

(¿Podriais ponerlo si sabeis en formato md? Hay simbolos que sé poner y otros no...)

Sea [a,b] un intervalo y P = {X0=a < X1 <...< Xn-1 < Xn=b} una partición del propio intervalo tal que {X0,X1...Xn} son los nodos de P
Se define el espacio de spline cuadráticos como el conjunto de funciones a trozos o splines definidas en el intervalo dado y asociadas a la partición P anterior tal que sus trozos son polinomios de grado menor o igual que dos de la forma ax^2 + bx + c, además, son funciones continuas y derivables en [a,b] (con derivada continua), es decir, son de clase 1, lo que nos proporcionará condiciones interesantes para resolver problemas de interpolantes con este tipo de splines.

Lo denotaremos como:  $S_2(x_1,x_2,...,x_n)$, tal que $(x_0,...,x_n)$ son los nodos de la partición P. Describamos a continuación este espacio.

En cuanto a su dimensión, es finita y esta es n+2. Se demuestra fácilmente:

- Cada trozo de un spline cuadrático (llamémoslo s) es de la forma ax^2 + bx + c, por lo tanto cada trozo está determinado por 3 parámetros. Supongamos que tenemos n trozos, por lo que tenemos 3n parámetros en total.

- Veamos las consecuencias de la continuidad y derivabilidad en todo el intervalo. Si imponemos esas condiciones tenemos que: $s_i(x_i)=s_{i+1}(x_{i+1})$ para todo i=1...n-1 y $s_i'(x_i)=s_{i+1}'(x_{i+1})$ para todo $i=1...n-1$. De cada condición se obtienen n-1 ecuaciones, por lo tanto obtendremos: n-1 + n-1 = 2n-2 ecuaciones linealmente independientes.

- Por lo tanto, $\dim(S_2(x_1,x_2,...,x_n))$ = 3n-(2n-2) = n+2, como queríamos demostrar.

Por otra parte, y con el conocimiento de la dimensión del espacio, podemos describir una base representativa del espacio de splines cuadráticos con el uso de potencias truncadas.

Una base del espacio es : $(1, x, x^2, (x-x_1)+^2, ... , (x-x_{n-1})+^2} con n+2 vectores linealmente independientes.


///// Lo de abajo es de Pablo

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


/* Conclusión de cuadráticos
   El problema de este procedimiento se presenta cuando hay que especificar condicinoes
   respecto a la derivada del interpolante en los puntos extremos x_0 y x_n. No
   existe un número suficiente de constantes para asegurar que se satisgafan las condiciones.
*/

###Splines cúbicos a partir de las segundas derivadas:
Uno de los problemas de la interpolación polinomial es que, al ir aumentando los
nodos (diferentes), el grado del polinomio aumenta (gr(p(x)) = nodos - 1). Esto
conlleva unas fluctuaciones en los extremos de la interpolación. (**1)

Sin embargo, si dividimos el intervalo en una partición P(n)={t_0=x_0<x_1<...<x_n=t_n}
, con un serie de subintervalos,
podemos aproximar un polinomio en cada intervalo minimizando la cota de error. (**3)

Para no
volver a tener el problema de las fluctuaciones, indeseables en la mayoría de las
aplicaciones, se suelen utilizar polinomios interpolantes de grado <= 3.
Esta técnica se conoce como aproximación polinomial fragmentaria, donde: (#1)

/* Dato curioso
   La palabra spline con el tiempo se usó para referirse a una larga banda flexible
   generalmente de metal, que podía usarse para dibujar curvas continuas suaves,
   forzando a la banda a pasar por puntos específicos y trazados a lo largo de la curva.
   (**2)
*/

// Este tipo de aproximación tiene una desventaja, no se tiene asegurada la derivabilidad

#### Spline cúbicos:
La aproximación más utilizada es la interpolación de splines cúbicos debido a que
proporciona un excelente ajuste a los puntos tabulados y su cálculo no es excesivamente
complejo.
Definimos una potencia truncada como: (**4)
Una potencia truncada pertenece a Clase k-1, su derivada k-1 es continua.

// Propiedades
Dada una función definida en [a,b], una partición del intervalo P(n)={a=x_0<x_1<...<x_n=b}
-S(x) es un polinomio cúbico denotado por S_j(x) en el subintervalo de extremos
x_j y x_(j+1), para j=0,1,..n-1.
-S_j(x_j) = f(x_j) y S_j(x_(j+1)) = f(x_(j+1))
-S'_(j+1)(x_(j+1)) = S'_j(x_(j+1))
-S''_(j+1)(x_(j+1)) = S''_j(x_(j+1))


Dentro de los cúbicos encontramos los de clase 1 y 2, denotados por: (#)
1.Los splines cúbicos de clase 1 son continuos y derivables en su dominio. Son un
espacio vectorial de dimensión 2*(n+1), cuya base es: (**5).

Una desventaja de estos splines es que no se asegura que haya derivabilidad en los
extremos, en un contexto geométrico eso significa que la función no es "suave" en
los puntos de unión. Generalmente las condiciones físicas necesitan esa suavidad,
y es aquí donde intervienen los splines cúbicos de clase 2.

2. Los splines cúbicos de clase 2 son continuos y 2 veces derivables.
Como sabemos que la dimensión de un spline S_k^r es (k-r)+r+1 la dimensión de 
este espacio vectorial es (3-2)n+2+1=n+3. Cuando k=r+1 el superíndice se omite.

Como tenemos n+1 variables, tenemos 2 libertades en la resolución.
Un tipo de spline es el Not-a-knot, requiere que la tercera derivada en los puntos
x_1 y x_(n-1) sea continua. Esto es S'''_0(x_1)=S'''_1(x_1) y 
S'''_(n-1)(x_(n-1))=S'''_n(x_(n-1)).





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





