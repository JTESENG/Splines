---
title: Ampliación de interpolación con Splines
author: [Miguel Anguita Ruiz, Pablo Baeyens Fernández, Pablo David Medina Sánchez, Rubén Morales Pérez, Francisco Javier Morales Piqueras]
lang: spanish
header-includes:
	\usepackage{mathrsfs}
	\usepackage{amsthm}
	\usepackage{booktabs}
	\usepackage{caption}
	\newtheorem*{proposicion}{Proposición}
	\newtheorem*{teorema}{Teorema}
	\theoremstyle{definition}
	\newtheorem*{definicion}{Definición}
	\newtheorem*{problema}{Problema}
	\theoremstyle{remark}
	\newtheorem*{solucion}{Solución}
toc: true
numbersections: true
fontsize: 11pt
geometry: margin=1in
---

\pagebreak

# Splines cuadráticos

##Introducción a los splines

\begin{definicion}
Sea $[a,b]$ un intervalo, $P = \{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$,
$k,r \in \mathbb{N}$, $r < k$ <!-->-->. Se dice que $s:[a,b] \to \mathbb{R}$ es un
spline si $s \in C^r([a,b])$ y para todo $1 \leq i \leq n$,
$s_{|[x_{i-1},x_i]} \in \mathbb{P}_k$. $S^r_k(P)$ es el espacio de dichas funciones.
\end{definicion}

La palabra **spline** con el tiempo se usó para referirse a una larga banda flexible
generalmente de metal, que podía usarse para dibujar curvas continuas suaves,
forzando a la banda a pasar por puntos específicos y trazados a lo largo de la curva.

## Descripción del espacio de splines cuadráticos

Partimos de $[a,b]$ un intervalo y $P \in \mathscr{P}([a,b])$. En esta primera sección
nos centramos en los splines cuadráticos: los pertenecientes a $S_2^1(P)$.

Sus trozos son polinomios de grado menor o igual que $2$ de la forma
$ax^2 + bx + c$. Además son funciones de clase $1$ (derivables en $[a,b]$
con derivada continua), lo que proporciona condiciones interesantes
para resolver problemas de interpolantes.

Veamos algunas propiedades:

\begin{proposicion}
Sea $[a,b]$ intervalo, $P = \{x_i\}_{i=0...n} \in \mathscr{P}([a,b])$, entonces $dim(S_2(P)) = n+2$.
\end{proposicion}

\begin{proof}

Sea $s \in S_2(P)$.

- Para cada intervalo $[x_{i-1}, x_i]$ $ s|_{[x_{i-1}, x_i]}(x) = ax^2 + bx + c$ para ciertos $a,b,c \in \mathbb{R}$. Por lo tanto cada trozo está determinado por 3 parámetros. Con $n$ trozos tenemos $3n$ parámetros en total.

- Si imponemos la continuidad y derivabilidad en los extremos tenemos que

\[
s_i(x_i)=s_{i+1}(x_i) \qquad s_i'(x_i)=s_{i+1}'(x_i)
\]

para todo $i=1...n-1$. De cada condición se obtienen $n-1$ ecuaciones, por
lo tanto obtendremos: $n-1 + n-1 = 2n-2$ ecuaciones linealmente
independientes.

Por lo tanto, $dim(S_2(P)) = 3n-(2n-2) = n+2$.
\end{proof}

Con el conocimiento de la dimensión del espacio podemos describir
una base del espacio de splines cuadráticos con el uso de potencias truncadas.
 Una **base del espacio** es: 

$$\{1, x, x^2, (x-x_1)_+^2, ... , (x-x_{n-1})_+^2\}$$.

\pagebreak

## Interpolación con splines cuadráticos

### Método local para la resolución de splines cuadráticos

El problema que debemos resolver es el siguiente:

\begin{problema}
Sea $[a,b]$ intervalo, $P \in \mathscr{P}([a,b])$ partición. Hallar $s \in S_2(P)$ tal que:
$$s(x_i)=y_i\ i=0,1,...,n$$
$$s'(x_k)=d_k$$
\end{problema}

Es decir, sabemos los valores de la función en todos los nodos, y sabemos el valor de la derivada en el nodo k.

\begin{solucion}
Vamos calcular para un nodo $i$ el polinomio $s_i(x)$, teniendo la condición $i>0$, y que sabemos la derivada en el punto: $d_i$.

\begin{table}[h]
\centering
\begin{tabular}{llll}
\hline
x & y & DD1 & DD2\\
\hline
$x_{i-1}$ & $y_{i-1}$ & & \\
$x_i$ & $y_i$ & $p_i=\frac{y_i-y_{i-1}}{h_i}$ & \\
$x_i$ & $y_i$ & $d_i$ & $\frac{d_i-p_i}{h_i}$\\
\hline
\end{tabular}
\end{table}

Siendo $h_i=x_i-x_{i-1}$.
De esta forma, para cada $s_i$ ya tendríamos una fórmula:

$$s_i(x)=y_{i-1}+p_k(x-x_{i-1})+\frac{d_k-p_k}{h_i}(x-x_{i-1})(x-x_k)$$

En el caso de que $i<>n$, y que sabiendo la derivada en el punto $d_i$, calculamos la tabla de Diferencias divididas así:

\begin{table}[h]
\centering
\begin{tabular}{llll}
\hline
x & y & DD1 & DD2\\
\hline
$x_i$ & $y_i$ & & \\
$x_i$ & $y_i$ & $d_i$ & \\
$x_{i+1}$ & $y_{i+1}$ & $p_{i+1}=\frac{y_{i+1}-y_i}{h_{i+1}}$ & $\frac{p_{i+1}-d_i}{h_i}$ \\
\hline
\end{tabular}
\end{table}

Siendo $s_i$ de la siguiente forma:

$s_i(x)=y_i+d_i(x-x_i)+\frac{p_{i+1}-d_i}{h_{i+1}}(x-x_i)(x-x_{i+1})$

Entonces, cuando empezamos, solo tenemos el nodo $k$, del que sabemos su derivada. De esta forma, usamos el método que haga falta.

Ahora, tendremos que que ir despejando el polinomio en las ecuaciones que que que no quedan de esta forma.

\begin{itemize}

\item para $i$ desde $k-1$ hasta $0$, decrementando $i$ en cada paso:
conocemos $s_i(x)$ así que evaluamos la derivadas en el nodo anterior
: $s'_{i}(x_{i-1})=d_{i-1}$, por lo tanto, tenemos $y_{i-1}$, $y_i$, y
 $d_{i-1}$, así que aplicamos el método primero para hallar $s_{i-1}(x)$.

\item para $i$ desde $k+1$ gasta $n$, incrementando $i$ en cada pasos:
sabemos que $s'_i(x_{i+1})=d_{i+1}$, tenemos $y_{i+1}$ $y_i$, y $d_{i+1}$, así que aplicamos el método segundo para hallar $s_{i+1}(x)$.
\end{itemize}

\end{solucion}

### Método global: cálculo con una base de potencias truncadas

Para este método usaremos una base del espacio vectorial $S_2(x_0,x_1...,x_n)$.

Tenemos los siguientes matrices y vectores:
$G$:matriz de Gram de nuestra base, en la cual evaluaríamos los elementos de la base en todos los nodos
$X$:vector de coeficientes
$b$:vector con los valores que queremos interpolar.

De esta forma, deberíamos resolver el sistema $G\ x=b$.


## Error en los splines cuadráticos

\begin{teorema}
Sean $f \in C^2([a,b])$, $\{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$,
$s \in S_2^1(\{x_i\}_{i = 0...n})$ spline para $f$,

$h = max\{x_i - x_{i-1}\}_{i = 1...n}$, $E = f - s$.
Además, sea $M >0$ tal que:

\[
 M \geq Sup \{|f''(x) - f''(y)| \; : \; |x - y| \leq h, \; x,y \in [a,b] \}
 \]

Entonces, se verifica, para todo $x \in [a,b]$:

\begin{equation}
E(x) \leq \frac{h^2M}{2}
\end{equation}

\end{teorema}

La demostración, así como cotas para las derivadas y cotas más precisas en
función de la localización de $x$ puede encontrarse en
*Quadratic Interpolatory Splines*,
W. Kammerer, G. Reddien y R.S. Varga, (1973).

## Ejemplos

\begin{problema}
Dados los datos de la tabla, halla mediante el método global el spline
cuadrático que interpole los nodos y cuya derivada en $x_1$ sea $4$.
\begin{table}[h]
\centering
\begin{tabular}{@{}l|lllll@{}}
$x_i$ & 2 & 4 & 5 & 8 \\
$y_i$ & 7 & 3 & 5 & 5 \\
$d_i$ &   & 4 &   &   \\
\end{tabular}
\end{table}
\end{problema}

\begin{solucion}

Debemos hallar $s \in S_2(P)$ con $a,b,c, \alpha, \beta \in \mathbb{R}$ tales que, para $x \in [2,8]$:

\[ s(x)  = a + bx + cx^2 + \alpha (x - 4)^2_+ + \beta (x -5)^2_+\]

Planteamos el sistema de ecuaciones $GX = b$:

\[
\begin{pmatrix}
1 & 2 & 4  & 0  & 0\\
1 & 4 & 16 & 0  & 0\\
1 & 5 & 25 & 0  & 0\\
1 & 8 & 64 & 16 & 9\\
0 & 1 &  8 & 0  & 0
\end{pmatrix}
\begin{pmatrix}
a 		 \\
b 		 \\
c 		 \\
\alpha \\
\beta
\end{pmatrix}
=
\begin{pmatrix}
7\\
3\\
5\\
5\\
4
\end{pmatrix}
\]

Resolviendo el sistema, obtenemos la solución $a = 35, b = -20, c = 3, \alpha = -5, \beta = 2$. Por tanto, para $x \in [2,8]$:

\[ s(x)  = 35 - 20x + 3x^2 -5(x - 4)^2_+ + 2(x -5)^2_+\]

Es decir:

\[
s(x) =
\begin{cases}
 3x^2 - 20x + 35  & x \in [2,4) \\
-2x^2 + 20x - 45 & x \in [4,5) \\
 5  & x \in [5,8] \\
\end{cases}
\]

\end{solucion}
\pagebreak

# Splines cúbicos

Uno de los problemas de la interpolación polinomial es que, al ir aumentando el
número de nodos el grado del polinomio necesario para interpolarlos aumenta.
Esto conlleva fluctuaciones en los extremos de la interpolación. <!--(1)-->

Si dividimos el intervalo en una partición podemos interpolar utilizando un
polinomio en cada intervalo, es decir, utilizando **splines cúbicos**. Como veremos después este método minimiza la cota de error.

\begin{equation}
	S(x) =
	\begin{cases}
	S_0(x) 			& \text{si } x \in {[t_0,t_1)} \\
	S_1(x)			& \text{si } x \in {[t_1,t_2)} \\
	\vdots			& \vdots \\
	S_{n-1}(x)		& \text{si } x \in {[t_{n-1},t_n)}
	\end{cases}
\end{equation}
 
Esta interpolación lineal fragmentaria pasa por los puntos: 
${ \{ (x_0,f(x_0)),(x_1,f(x_1)),...,(x_n,f(x_n)) \} }$


Dentro de los cúbicos encontramos los de clase 1 y 2, denotados por $S^{1}_3$ y $S^{2}_3$ (ó $S_3$).

1. Los splines cúbicos de clase 1 son continuos y derivables
con derivada continua. Forman un espacio vectorial de dimensión $2(n+1)$, cuya base es:

${ \{1,x,x^2,x^3, (x-x_1)^{2}_{+},(x-x_1)^{3}_{+},...,(x-x_{n-1})^{2}_{+},(x-x_{n-1})^3_+ \} }$

Estos splines no aseguran derivabilidad en los extremos.
En un contexto geométrico esto significa que la función no es *suave* en
los puntos de unión. Generalmente las condiciones físicas necesitan esa suavidad,
y es aquí donde intervienen los splines cúbicos de clase 2.

2. Los splines cúbicos de clase 2 son continuos y 2 veces derivables.
Como sabemos que la dimensión de un spline la dimensión de este espacio es
$(3-2)n+2+1=n+3$.

Como tenemos $n+1$ variables, tenemos $2$ libertades en la resolución.

## Construcción a partir de los valores de $s''$ en los nodos $\{x_i\}$

Vamos a plantear un método de resolución utilizando las segundas derivadas, denotamos,
para $i=1, ... n-1$:

- $M_i = S''(x_i)$, que son desconocidos a priori salvo en un spline natural.
- $h_i = x_i-x_{i-1}$

Como el spline es de clase 2, tenemos para $i=1, ... {n-1}$:
$$S''(x_i) = S''_i(x_i) = S''_{i+1}(x_i)$$

La restricción a cada intervalo de $S$ es un polinomio $S_i$ de grado 3, por ende, $S''_i$ es lineal, con expresión para $x \in [x_{i-1},x_i]$:


$$S''_i(x) = M_{i-1} \frac{x_i-x}{h_i} + M_i\frac{x-x_{i-1}}{h_i}$$

Integramos dos veces usando que $S_i(x_{i-1}) = y_{i-1}$ y $S_i(x_i) = y_i$ para las constantes de integración obteniendo, para ${x \in [x_{i-1},x_i]}$:

$$S_i(x) = M_{i-1}\frac{(x_i-x)^3}{6h_i} + M_i\frac{x-x_{i-1}}{6h_i} + \frac{y_{i-1} - M_{i-1}h^2_i}{6} \cdot \frac{x_i-x}{h_i} + \frac{y_i- M_ih^2_i}{6} \cdot \frac{x-x_{i-1}}{h_i}$$ 

Esta ecuación nos permite calcular $S$ conocidas $M_i$ con $i=0,1,...n$.
Las condiciones de suavidad en las ligaduras nos permiten igualar $S'_{i+1}(x_i) = S'_i(x_i)$. Derivando una vez, si $x \in {[x_{i-1},x_i]}$:

$$S'_i(x) = -M_{i-1}\frac{(x_i-x)^2}{2h_i} + M_i\frac{(x-x_{i-1})^2}{2h_i} + \frac{y_i-y_{i-1}}{h_i} -(M_i-M_{i-1})\frac{h_i}{6}$$

Si $x \in {[x_{i},x_{i+1}]}$: 

$$S'_{i+1}(x) = -M_i\frac{(x_{i+1}-x)^2}{2h_i} + M_{i+1}\frac{(x-x_i)^2}{2h_{i+1}}
 + \frac{y_{i+1}-y_i}{h_{i+1}} -(M_{i+1}-M_i)\frac{h_{i+1}}{6}$$


Recordando que $h_i = x_i - x_{i-1}$ e igualando $S'_{i+1}(x_i) = S'_i(x_i)$:

$$-M_i\frac{h_{i+1}}{2} + \frac{y_{i+1}-y_i}{h_{i+1}} -(M_{i+1}-M_i)\frac{h_{i+1}}{6}
=  M_i\frac{h_i}{2} + \frac{y_i-y_{i-1}}{h_i} -(M_i-M_{i-1})\frac{h_i}{6}$$ 

Agrupamos los $M_i$: 

$$-M_i\frac{h_{i+1}}{2} + M_i\frac{h_{i+1}}{6} - M_i\frac{h_i}{2} + M_i\frac{h_i}{6} + \frac{y_{i+1}-y_i}{h_{i+1}} - \frac{y_i-y_{i-1}}{h_i} =  M_{i+1}\frac{h_{i+1}}{6} + M_{i-1}\frac{h_i}{6}$$

Multiplicamos a ambos lados por $6$, sacamos factor común y recordamos que $f[x_i,x_{i+1}] = \frac{y_{i+1}-y_i}{h_{i+1}}$:

$$6M_i\frac{-3h_{i+1}}{6} + \frac{h_{i+1}}{6} - 3\frac{h_i}{6} + \frac{h_i}{6} + 6(f[x_i,x_{i+1}] - f[x_{i-1},x_i]) =  M_{i+1}h_{i+1} + M_{i-1}h_i$$


Agrupando y multiplicando $M_i$ arriba y abajo por $-2$: 

$$-2M_i\frac{-2h_{i+1}-3h_i+h_i}{-2} + 6(f{[x_i,x_{i+1}]} - f{[x_{i-1},x_i]}) =  M_{i+1}h_{i+1} + M_{i-1}h_i$$

Pasamos el $M_i$ a la derecha y dividimos por $(h_{i+1}+h_i)$ en ambos lados:

$$6\frac{f{[x_i,x_{i+1}]} - f{[x_{i-1},x_i]}}{h_{i+1}-h_i} =  M_{i+1}\frac{h_{i+1}}{h_{i+1}+h_i} + M_{i-1}\frac{h_i}{h_{i+1}+h_i} + 2M_i$$

$$6f{[x_{i-1},x_i,x_{i+1}]} = \frac{M_{i+1}h_{i+1} + M_{i-1}h_i}{h_{i+1}+h_i} + 2M_i$$


Denotando por $m_i = \frac{h_i}{h_i+h_{i+1}}$, $\lambda_i = \frac{h_{i+1}}{h_i+h_{i+1}}$ y $\gamma_i = 6f[x_{i-1},x_i,x_{i+1}]$:

$$m_iM_{i-1} + 2M_i + \lambda_iM_{i+1} = \gamma_i$$

Con los $M_i$ en las ligaduras tendremos $4(n-1)$ variables, para que el sistema
sea determinado nos faltan dos condiciones. Hay diferentes condiciones que se nos pueden presentar:

\vspace*{2\baselineskip}

**Spline sujeto**
 
$S'_1(x_0) = f'_0$ y $S'_n(x_n)=f'_n$. De acuerdo con la fórmula de $S'(x)$ obtenemos:

$$f'_0 = -\frac{M_0h_i}{2} + f[x_0,x_1] - \frac{(M_1 - M_0)h_i}{6} \implies  2M_0+M_1=\frac{6(f{[x_0,x_1]} - f^{'}_0)}{h_1} = 6f{[x_0,x_0,x_1]}$$

Equivalentemente para $x_n$:

\begin{multline*}
S'_n(x_n) = - \frac{M_{n-1}(x_n-x_n)^2}{2h_n} + \frac{M_n(x_n-x_{n-1})^2}{2h_n} + \frac{(y_n-y_{n-1})}{h_n} - \frac{(M_n-M_{n-1})h_n}{6} \\ \implies
M_{n-1}+2M_n=6f[x_{n-1},x_n,x_n]
\end{multline*}

Tomando $\mu_i = h_i/(h_i+h_{i+1})$, la matriz del sistema es:


$$
\begin{pmatrix}
  2 	   & \lambda_0 &    0       &   \cdots  &     0	         \\
  \mu_1  & 2	 		& \lambda_1  &   0       &    \vdots      \\
  0      & \ddots    & \ddots     &  \ddots   &     0          \\
  \vdots &     0     & \mu_{n-1}  &    2      & \lambda_{n-1}  \\
  0      &   \cdots  &     0      &   \mu_n   &     2
\end{pmatrix}
\begin{pmatrix}
  M_0 \\
  M_1 \\
  \vdots \\
  M_{n-1} \\
  M_n
\end{pmatrix} =
\begin{pmatrix}
  \gamma_0 \\
  \gamma_1 \\
  \vdots \\
  \gamma_{n-1} \\
  \gamma_n
\end{pmatrix}$$

\vspace*{2\baselineskip}

**Spline natural** 

En este caso $M_0=0$ y $M_n=0$, por lo que el sistema queda:

$$\begin{pmatrix}
  2 	   & \lambda_0 &    0       &   \cdots  &     0	         \\
  \mu_1  & 2	 		& \lambda_1  &   0       &    \vdots      \\
  0      & \ddots    & \ddots     &  \ddots   &     0          \\
  \vdots &     0     & \mu_{n-1}  &    2      & \lambda_{n-1}  \\
  0      &   \cdots  &     0      &   \mu_n   &     2
\end{pmatrix}
\begin{pmatrix}
  M_0 \\
  M_1 \\
  \vdots \\
  M_{n-1} \\
  M_n
\end{pmatrix}
\begin{pmatrix}
  0 \\
  \gamma_1 \\
  \vdots \\
  \gamma_{n-1} \\
  0
\end{pmatrix}$$

\vspace*{2\baselineskip}

**Spline periódico**

En este caso $S'_1(x_0) = S'_n(x_n)$ y $S''_1(x_0) = S''_n(x_n)$. El sistema queda:

<!--Falta plantearlo-->

$$
\begin{pmatrix}
  2 	   & \lambda_0 &    0       &   \cdots  &     0	         \\
  \mu_1  & 2	 		& \lambda_1  &   0       &    \vdots      \\
  0      & \ddots    & \ddots     &  \ddots   &     0          \\
  \vdots &     0     & \mu_{n-1}  &    2      & \lambda_{n-1}  \\
  0      &   \cdots  &     0      &   \mu_n   &     2
\end{pmatrix}
\begin{pmatrix}
  M_0 \\
  M_1 \\
  \vdots \\
  M_{n-1} \\
  M_n
\end{pmatrix} =
\begin{pmatrix}
  0 \\
  \gamma_1 \\
  \vdots \\
  \gamma_{n-1} \\
  0
\end{pmatrix}$$

**CAMBIAR SISTEMA POR EL QUE ES!!!**

## Propiedades de minimización

Comenzamos planteando un problema de minimización sobre el espacio euclídeo
 $(C^2([a,b]), <\cdot,\cdot>)$, con la métrica y norma definida de la forma usual:

$$<f,g> = \int_a^b fg, \qquad || f || = \sqrt{ \int_a^b f^2}$$

\vspace*{\baselineskip}

Planteamos el problema:

\begin{problema}
Sea $f \in C^2([a,b])$, $P \in \mathscr{P}([a,b])$. Sea $H \subset C^2([a,b])$ definido por:
\[H = \{g \in C^2([a,b]) \; : \; \forall p \in P \; g(p) = f(p) \text{ y } \; g'(a) = f'(a), \; g'(b) = f'(b)\} \]

Hallar $u \in H$ tal que $||u''||$ sea mínima.
\end{problema}

\vspace*{\baselineskip}

Para resolver el problema, demostramos el siguiente teorema:

\begin{teorema}[Minimización]
Sea $f \in C^2([a,b])$, $P \in \mathscr{P}([a,b])$, $s$ spline sujeto para $f$. Se verifica:

\[
\forall u \in H : \; ||s''|| \leq ||u''||
\]
\end{teorema}

\vspace*{\baselineskip}

\begin{proof}
Sea $u \in H$, $e = u - s$. Tenemos:

\[
||u''||^2 = ||e'' + s''||^2 = ||e''||^2 + ||s''||^2 + 2<e'', s''>
\]

Dividimos $<e'',s''>$ en intervalos:
\[
<e'',s''> = \int_a^b e''s''
= \sum_1^{n-1} \int_{x_i}^{x_{i+1}} e''s''
\]

En cada intervalo, integramos por partes:
\[
\sum_1^{n-1}  \int_{x_i}^{x_{i+1}} e''s''
= \sum_1^{n-1}  \left. e'(x)s''(x) \right|_{x_i}^{x_{i+1}} - \sum_1^{n-1}  \int_{x_i}^{x_{i+1}}e's'''
\]

La primera sumatoria es una suma telescópica, por lo que conservamos
el primer y último término:
\[
\sum_1^{n-1} \left. e'(x)s''(x) \right|_{x_i}^{x_{i+1}}
= e'(b)s''(b) - e'(a)s''(a)
= (u'(b) - s'(b))s''(b) - (u'(a) - s'(b))s''(a)
= 0
\]

ya que $u, s \in H$.

En cuanto a la segunda, $s'''|_{[x_i, x_{i+1}]}$ es constante, por lo
que podemos sacarlo de la integral:
\[
\sum_{1}^{n-1} s_i\int_a^be'(x) = \sum_{1}^{n-1}  s_i(e(b) - e(a)) = 0
\]

Es decir, $<e'',s''> = 0$. Por tanto:

\[
||u''||^2 = ||e''||^2 + ||s''||^2 + 2<e'', s''>
= ||e''||^2 + ||s''||^2 \geq ||s''||^2
\]

donde utilizamos que la norma siempre es positiva.
\end{proof}

Así, podemos observar que el **spline cúbico sujeto**
asociado a una función $f$ tiene la menor norma de su segunda derivada
de entre las que interpolan a $f$ en una partición dada, por lo que
resuelve nuestro problema.

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
