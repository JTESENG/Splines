---
title: Ampliación de interpolación con Splines
author: [Miguel Anguita Ruiz, Pablo Baeyens Fernández, Pablo David Medina Sánchez, Rubén Morales Pérez, Francisco Javier Morales Piqueras]
lang: spanish
header-includes:
	\usepackage{mathrsfs}
	\usepackage{amsthm}
	\usepackage{booktabs}
	\usepackage{caption}
	\usepackage{xfrac}
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

La palabra **spline** con el tiempo se usó para referirse a una larga banda flexible
generalmente de metal, que podía usarse para dibujar curvas continuas suaves,
forzando a la banda a pasar por puntos específicos y trazados a lo largo de dicha curva.

La formalización del concepto de función spline, es decir, una curva continua
que pasa por ciertos puntos se resume en la siguiente definición:

\begin{definicion}
Sea $[a,b]$ un intervalo, $P = \{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$,
$k,r \in \mathbb{N}$, $r < k$. Se dice que $s:[a,b] \to \mathbb{R}$ es un
spline si $s \in C^r([a,b])$ y para todo $1 \leq i \leq n$,
$s_{|[x_{i-1},x_i]} \in \mathbb{P}_k$. $S^r_k(P)$ es el espacio de dichas funciones.
\end{definicion}

## Descripción del espacio de splines cuadráticos

Partimos de $[a,b]$ un intervalo y $P \in \mathscr{P}([a,b])$. En esta primera sección
nos centramos en los splines cuadráticos: los pertenecientes a $S_2^1(P)$.

Sus trozos son polinomios de grado menor o igual que $2$ de la forma
$ax^2 + bx + c$. Además son funciones de clase $1$ (derivables en $[a,b]$
con derivada continua), lo que proporciona unas condiciones interesantes
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

$$\{1, x, x^2, (x-x_1)_+^2, ... , (x-x_{n-1})_+^2\}$$

\pagebreak

## Interpolación con splines cuadráticos

### Método local: cálculo trozo a trozo

El problema que debemos resolver es el siguiente:

\begin{problema}
Sea $[a,b]$ intervalo, $P \in \mathscr{P}([a,b])$ partición. Hallar $s \in S_2(P)$ tal que:
$$s(x_i)=y_i\ i=0,1,...,n$$
$$s'(x_k)=d_k$$
\end{problema}

Es decir, sabemos los valores de la función en todos los nodos y el valor de la derivada en el nodo $k$.

\vspace*{2\baselineskip}

\begin{solucion}
Si $k > 0$, para calcular $s_k$ podemos calcular la tabla de diferencias divididas:

\vspace*{2\baselineskip}

\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$x_{k-1}$ & $y_{k-1}$ & & \\
$x_k$ & $y_k$ & $p_k$ & \\
$x_k$ & $y_k$ & $d_k$ & $\frac{d_k-p_k}{h_k}$\\
\hline
\end{tabular}
\end{table}

\vspace*{2\baselineskip}

De esta forma, $s_k$ queda, para $x \in [x_{k-1}, x_k]$

\begin{equation} \label{eq:sk}
s_k(x)=y_{k-1}+p_k(x-x_{k-1})+\frac{d_k-p_k}{h_k}(x-x_{k-1})(x-x_k)
\end{equation}

Conocida la expresión de $s_k$ podemos calcular $d_{k-1} = s'_k(x_{k-1})$, y repetir
este proceso para calcular $s_{k-1}$, hasta llegar a $k = 0$.


Si $k < n$, debemos calcular $s_{k+1}$. Como sabemos la derivada $d_k$, calculamos la tabla de diferencias divididas:

\vspace*{2\baselineskip}

\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$x_k$ & $y_k$ & & \\
$x_k$ & $y_k$ & $d_k$ & \\
$x_{k+1}$ & $y_{k+1}$ & $p_{k+1}$ & $\frac{p_{k+1}-d_k}{h_k}$ \\
\hline
\end{tabular}
\end{table}

\vspace*{2\baselineskip}

De esta forma, $s_{k+1}$ queda para $x \in [x_k, x_{k+1}]$ de la siguiente forma:

\begin{equation} \label{eq:skmas}
s_{k+1}(x)=y_k+d_k(x-x_k)+\frac{p_{k+1}-d_k}{h_{k+1}}(x-x_k)(x-x_k)
\end{equation}


El método queda entonces de la siguiente forma:

\begin{enumerate}

\item Para $i$ desde $k$ hasta $0$:
\begin{enumerate}
\item Calculamos $d_i$, (conocida en el primer caso) haciendo $d_i = s'_{i+1}(x_i)$.
\item Aplicamos la fórmula (\ref{eq:sk}) para calcular $s_i$.
\end{enumerate}

\item Para $i$ desde $k+1$ hasta $n$:
\begin{enumerate}
\item Calculamos $d_i$ haciendo $d_i = s'_{i-1}(x_{i})$
\item Aplicamos la fórmula (\ref{eq:skmas}) para calcular $s_i$.
\end{enumerate}
\end{enumerate}

\end{solucion}

\vspace*{2\baselineskip}

### Método global: cálculo con una base de potencias truncadas

Para este método usaremos esta base del espacio vectorial $S_2(P)$:

$$\{1, x, x^2, (x-x_1)_+^2, ... , (x-x_{n-1})_+^2\}$$

Tenemos los siguientes matrices y vectores:

- $G$: matriz de Gram. Evaluamos los elementos de la base en todos los nodos.
- $X$: vector de coeficientes
- $b$: vector con los valores que queremos interpolar.

De esta forma, deberíamos resolver el sistema $GX=b$.

Si notamos por $x_k$ el nodo en el que conocemos la derivada y $d_k$ la derivada
en el nodo, el sistema queda:

\begin{equation*}
\begin{pmatrix}
	1 & x_0 & x_0^2   & 0 & \cdots & 0\\
	1 & x_1 & x_1^2   & (x_1-x_1)_{+}^2 & \cdots & 0\\
	\vdots& & \vdots  & \vdots          & \cdots & \vdots \\
	\vdots& & \vdots  & \vdots          & \cdots & \vdots \\
	1 & x_n & x_n^2   & (x_n-x_1)_{+}^2 & \cdots & (x_n-x_{n-1})_{+}^2\\
	0 &   1 &  2x_k & 2(x_k-x_1)_{+} & \cdots & 2(x_k-x_{n-1})_{+}
\end{pmatrix}
\begin{pmatrix}
	a \\
	b \\
	c \\
	\alpha \\
	\vdots \\
	\omega
\end{pmatrix}
=
\begin{pmatrix}
	y_0\\
	y_1\\
	\vdots\\
	\vdots\\
	y_n\\
	d_k\\
\end{pmatrix}
\end{equation*}

Es una matriz escalonada ya que las potencias truncadas serán $0$ antes del nodo
que las define. Finalmente la solución sería, para $x \in [a,b]$:

$$s(x) = a + bx + cx^2 + \alpha(x-x_1)_+^2 + \cdots + \omega(x-x_{n-1})_+^2$$

## Error en los splines cuadráticos

\begin{teorema}
Sean $f \in C^2([a,b])$, $\{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$,
$s \in S_2^1(\{x_i\}_{i = 0,...,n})$ spline para $f$,

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
cuadrático que interpole los nodos $x_i$ con $i=0,...,3$ y cuya derivada en $x_1$ sea $4$.
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
 3x^2 - 20x + 35  & x \in [2,4] \\
-2x^2 + 20x - 45 & x \in (4,5] \\
 5  & x \in (5,8] \\
\end{cases}
\]

\end{solucion}

\begin{problema}
Dados los siguientes dados, calcula el spline cuadrático que los interpola:
\begin{table}[h]
\centering
\begin{tabular}{l|lllll}
$x_i$ & $-1$& $1$ & $3$ & $6$ & $7$ \\
$y_i$ & $1$ & $4$ & $8$ & $2$ & $9$ \\
$d_i$	&     &	   & $5$ &		&     \\
\end{tabular}
\end{table}
\end{problema}

\begin{solucion}
Nos dan la derivada en el nodo 3, procedemos a calcular las diferencias divididas en los nodos 1 y 3 para hallar $s_2$:


\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$1$ & $4$ &     & \\
$3$ & $8$ & $2$ & \\
$3$ & $8$ & $5$ & $\sfrac{3}{2}$ \\
\hline
\end{tabular}
\end{table}

$s_2$ queda en su intervalo:

$$s_2(x)=4+2(x-1)+\frac{3}{2}(x-1)(x-3)$$

Ahora estimamos la derivada en el nodo 1:

$$s_2'(x)=2+\frac{3}{2}((x-3)+(x-1))=3x - 4$$

$$s_2'(1)= 3 \cdot 1 - 4 = -1$$

Realizamos de nuevo la tabla de diferencias divididas:

\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$-1$ & $1$  &                & \\
$1$  & $4$  & $\sfrac{3}{2}$  & \\
$1$  & $4$  & $-1$           & $\sfrac{-5}{4}$ \\
\hline
\end{tabular}
\end{table}

$s_1$ queda en su intervalo:

$$s_1(x)=1+\frac{3}{2}(x+1)+\frac{-5}{4}(x+1)(x-1)$$

Ahora que hemos calculado la expresión de $s$ para todos lo intervalos a la izquierda de la derivada, calculamos la función para todos los valores a la derecha de la
derivada.

Calculamos las diferencias divididas para nodos 3 y 6:

\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$3$  & $8$  &      & \\
$3$  & $8$  & $5$  & \\
$6$  & $2$  & $-2$ & $\sfrac{-7}{3}$ \\
\hline
\end{tabular}
\end{table}




$s_3$ y su derivada quedan en su intervalo:
$$s_3(x)=8+5(x-3)- \frac{7}{3} (x-3)(x-3)$$

$$s_3'(x)=5-\frac{7}{3}2(x-3)=5-\frac{7}{3}(2x-9)$$

Estimamos la derivada del nodo 6:
$$s_3'(6)=5-\frac{7}{3}3 = -2$$

Finalmente, calculamos $s_4$:

\begin{table}[h]
\centering
\begin{tabular}{llll}
x & y & DD1 & DD2\\
\hline
$6$ & $2$  &       & \\
$6$ & $2$  & $-9$  & \\
$7$ & $9$  & $7$   & $16$ \\
\hline
\end{tabular}
\end{table}

De esta forma, la expresión de $s_4$ sería:

$$s_4(x)=2 -9(x-6)+16(x-6)(x-6)$$

Por lo tanto, nuestra solución sería:


$$s(x)=
\begin{cases}
s_1(x)=1+\frac{3}{2}(x+1)+\frac{-5}{4}(x+1)(x-1)   & \text{si } x\in {[-1,1)}\\
s_2(x)=4+2(x-1)+\frac{3}{2}(x-1)(x-3)             & \text{si } x\in {[1,3)}\\
s_3(x)=8+5(x-3)- \frac{7}{3} (x-3)(x-3)           & \text{si } x\in {[3,6)}\\
s_4(x)=2 -9(x-6)+16(x-6)(x-6) & \text{si } x\in {[6,7])} \\
\end{cases}
$$
\end{solucion}

\pagebreak

# Splines cúbicos

Uno de los problemas de la interpolación polinomial es que, al ir aumentando el
número de nodos el grado del polinomio requerido para interpolarlos aumenta.
Esto conlleva fluctuaciones en los extremos de la interpolación. <!--(1)-->

Si dividimos el intervalo en una partición podemos interpolar utilizando un
polinomio S_i(x) de grado 3 en cada intervalo, es decir, utilizando **splines cúbicos**. Como veremos después este método minimiza la cota de error.

\begin{equation}
	S(x) =
	\begin{cases}
	S_0(x) 			& \text{si } x \in {[x_0,x_1)} \\
	S_1(x)			& \text{si } x \in {[x_1,x_2)} \\
	S_i(x)   		& \text{si } x \in {[x_i,x_{i+1})} \\
	S_{n-1}(x)		& \text{si } x \in {[x_{n-1},x_n]}
	\end{cases}
\end{equation}

Esta interpolación lineal fragmentaria pasa por los puntos:
${ \{ (x_0,f(x_0)),(x_1,f(x_1)),...,(x_n,f(x_n)) \} }$


Dentro de los cúbicos encontramos los de clase 1 y 2, denotados por $S^{1}_3$ y $S^{2}_3$ (ó $S_3$).

1. Los splines cúbicos de **clase 1** son continuos y derivables
con derivada continua. Conforman un espacio vectorial de dimensión $2(n+1)$. Una base es:
$$\{1,x,x^2,x^3, (x-x_1)^{2}_{+},(x-x_1)^{3}_{+},...,(x-x_{n-1})^{2}_{+},(x-x_{n-1})^3_+ \}$$
Estos splines no aseguran derivabilidad en los extremos.
En un contexto geométrico esto significa que la función no es *suave* en
los puntos de unión. Generalmente las condiciones físicas necesitan esa suavidad,
y es aquí donde intervienen los splines cúbicos de clase 2.

2. Los splines cúbicos de **clase 2** son continuos y 2 veces derivables.
A partir de la fórmula general, la dimensión de este espacio para una partición
$\{x_i\}_{i=0,...,n}$ es $dim (S_3^2(P)) = (3-2)n+2+1=n+3$. Como tenemos $n+1$ variables,
tenemos $2$ libertades en la resolución.

## Construcción a partir de los valores de $s''$ en los nodos $\{x_i\}$

Vamos a plantear un método de resolución utilizando las segundas derivadas, denotamos,
para $i=1, ..., n-1$: $M_i = S''(x_i)$, que son desconocidos a priori salvo en un spline natural.

Como el spline es de clase 2, tenemos para $i=1, ... {n-1}$:
$$S''(x_i) = S''_i(x_i) = S''_{i+1}(x_i)$$

La restricción a cada intervalo de $S$ es un polinomio $S_i$ de grado 3, por ende, $S''_i$ es lineal, con expresión para $x \in [x_{i-1},x_i]$:


$$S''_i(x) = M_{i-1} \frac{x_i-x}{h_i} + M_i\frac{x-x_{i-1}}{h_i}$$

Integramos dos veces usando que $S_i(x_{i-1}) = y_{i-1}$ y $S_i(x_i) = y_i$ para las constantes de integración, obteniendo, para ${x \in [x_{i-1},x_i]}$:

$$S_i(x) = M_{i-1}\frac{(x_i-x)^3}{6h_i} + M_i\frac{x-x_{i-1}}{6h_i} + (y_{i-1}-\frac{ M_{i-1}h^2_i}{6}) \cdot \frac{x_i-x}{h_i} + (y_i-\frac{ M_ih^2_i}{6}) \cdot \frac{x-x_{i-1}}{h_i}$$

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

Multiplicamos a ambos lados por $6$, sacamos factor común y recordamos que $p_{i+1} = \frac{y_{i+1}-y_i}{h_{i+1}}$:

$$6M_i\frac{-3h_{i+1}}{6} + \frac{h_{i+1}}{6} - 3\frac{h_i}{6} + \frac{h_i}{6} + 6(p_{i+1} - p_i) =  M_{i+1}h_{i+1} + M_{i-1}h_i$$


Agrupando y multiplicando $M_i$ arriba y abajo por $-2$:

$$-2M_i\frac{-2h_{i+1}-3h_i+h_i}{-2} + 6(p_{i+1} - p_i) =  M_{i+1}h_{i+1} + M_{i-1}h_i$$

Pasamos el $M_i$ a la derecha y dividimos por $(h_{i+1}+h_i)$ en ambos lados:

$$6\frac{p_{i+1}-p_i}{h_{i+1}+h_i} =  M_{i+1}\frac{h_{i+1}}{h_{i+1}+h_i} + M_{i-1}\frac{h_i}{h_{i+1}+h_i} + 2M_i$$


Denotando por $\displaystyle\mu_i = \frac{h_i}{h_i+h_{i+1}}$, $\displaystyle\lambda_i = 1-\mu_i = \frac{h_{i+1}}{h_i+h_{i+1}}$ y $\displaystyle\gamma_i = 6\frac{p_{i+1}-p_i}{h_{i+1}+h_i}$:

$$\mu_iM_{i-1} + 2M_i + \lambda_iM_{i+1} = \gamma_i (*)$$

Con los $M_i$ en las ligaduras tendremos $4(n-1)$ variables, para que el sistema
sea determinado nos faltan dos condiciones. Hay diferentes condiciones que se nos pueden presentar:

\vspace*{2\baselineskip}

 - **Spline sujeto**

$S'_1(x_0) = f'_0$ y $S'_n(x_n)=f'_n$. De acuerdo con la fórmula de $S'(x)$ obtenemos:

$$f'_0 = -\frac{M_0h_i}{2} + f[x_0,x_1] - \frac{(M_1 - M_0)h_i}{6} $$

$$\implies  2M_0+M_1=\frac{6(f{[x_0,x_1]} - f^{'}_0)}{h_1} = 6f{[x_0,x_0,x_1]} (*)$$

Equivalentemente para $x_n$:

\begin{multline*}
S'_n(x_n) = - \frac{M_{n-1}(x_n-x_n)^2}{2h_n} + \frac{M_n(x_n-x_{n-1})^2}{2h_n} + \frac{(y_n-y_{n-1})}{h_n} - \frac{(M_n-M_{n-1})h_n}{6} \\
\implies
M_{n-1}+2M_n=6f[x_{n-1},x_n,x_n] (*)
\end{multline*}

Usando (*), la matriz del sistema es:


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

 - **Spline natural**

En este caso $M_0=0$ y $M_n=0$, $\lambda_0 = \mu_n = 1$ por lo que el sistema queda:

$$\begin{pmatrix}
   2        & \lambda_1  &   0       &    \cdots     & 0 \\
  \mu_2     & 2          &  \lambda_2   &     0         & \vdots\\
  0         & \ddots     &  \ddots   &     \ddots    & 0 \\
  \vdots    &    0       & \mu_{n-2} &     2         &  \lambda_{n-2}\\
  0         &  \cdots    &    0      &    \mu_{n-1}  & 2
\end{pmatrix}
\begin{pmatrix}
  M_1 \\
  M_2 \\
  \vdots \\
  M_{n-2} \\
  M_{n-1} \\
\end{pmatrix}
=
\begin{pmatrix}
  \gamma_1 \\
  \gamma_2 \\
  \vdots \\
  \gamma_{n-2} \\
  \gamma_{n-1} \\
\end{pmatrix}$$

\vspace*{2\baselineskip}

 - **Spline periódico**

En este caso $S'_1(x_0) = S'_n(x_n)$ y $S''_1(x_0) = S''_n(x_n)$. El sistema queda:

<!--Falta plantearlo-->

$$
\begin{pmatrix}
  2 	   & \lambda_0 &    0       &   \cdots  &     0	         \\
  \mu_1  & 2	 		& \lambda_1  &   0       &    \vdots      \\
  0      & \ddots    & \ddots     &  \ddots   &     0          \\
  \vdots &     0     & \mu_{n-1}  &    2      & \lambda_{n-1}
\end{pmatrix}
\begin{pmatrix}
  M_0 \\
  M_1 \\
  \vdots \\
  M_{n-2} \\
  M_{n-1}
\end{pmatrix} =
\begin{pmatrix}
  \gamma_0=h_1-S^{'}_(x_0) \\
  \gamma_1 \\
  \vdots \\
  \gamma_{n-2} \\
  \gamma_{n-1}
\end{pmatrix}$$

En este caso añadimos:

${ S^{'}_1(x_0)= -M_0 \cdot \frac{h_1}{2}+ f{[x_o,x_1]} - \frac{M_1-M_0}{6} \cdot h_1 }$
<!--**CAMBIAR SISTEMA POR EL QUE ES!!!**-->

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
\sum_{1}^{n-1} s_i\int_a^be' = \sum_{1}^{n-1}  s_i(e(b) - e(a)) = 0
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

## Error en los splines cúbicos

\begin{teorema}
Sea $f \in C^4([a,b])$, $n \in \mathbb{N}$, $P = \{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$ y  $s \in S_3^1(P)$ spline para $f$. Además, sean $h = max\{x_i - x_{i-1}\}_{i = 1...n}$, $M > 0$ cota superior de $|f^{iv)}|$ en $[a,b]$ y $E = f - s$, $x \in [a,b]$.

Se verifica:

\begin{equation} \label{eq:errorS31}
|E(x)| \leq \frac{5M}{384}h^4
\end{equation}

\end{teorema}

La demostración, así como cotas para las derivadas, puede consultarse en *Optimal Error Bounds for Cubic Spline Interpolation*, Charles Hall y Weston Meyer, (1976).

## Ejemplos

**Sujeto**:

\begin{problema}

Hallar spline sujeto tal que:
\begin{enumerate}
\item Pasa por los puntos: $\{(0,0),(1,0.5),(2,2), (3,1.5)\}$
\item $S'(0) = 0.2$ y $S'(3) = -1$
\end{enumerate}
\end{problema}

\vspace*{\baselineskip}

\begin{solucion}
Como los nodos están equiespaciados $h_i=1$ $\forall i \in \{1..n\}$

$$\lambda_0 = \mu_3 = 1 \; \text{ y } \; \lambda_1=\lambda_2=\mu_1=\mu_2= \frac{1}{2}$$

Calculamos las diferencias divididas para obtener los $\gamma_i$

\begin{itemize}

\item $\displaystyle\frac{\gamma_0}{6} = f[x_0,x_0,x_1] =  \frac{f[x_1,x_0]-f[x_0,x_0]}{x_1-x_0} = \left(\frac{0.5-0}{1-0}-0.2\right)/1 = 0.3$

\item $\displaystyle\frac{\gamma_1}{6} = f{[x_0,x_1,x_2]} =  \frac{f{[x_2,x_1]}-f{[x_1,x_0]} }{ x_2-x_0 } = \left(\frac{2-0.5}{1-0}-\frac{0.5-0}{1-0}\right)/2 = \frac{1}{2}$

\item $\displaystyle\frac{\gamma_2}{6} = f{[x_1,x_2,x_3]} =  \frac{ f{[x_3,x_2]}-f{[x_2,x_1]} }{ x_3-x_1 } = \left(\frac{1.5-2}{1-0}-\frac{2-0.5}{1-0}\right)/1 = -1$

\item $\displaystyle\frac{\gamma_3}{6} = f{[x_2,x_3,x_3]} =  \frac{ f{[x_3,x_3]}-f{[x_2,x_3]} }{ x_3-x_2 } = \left(-1-\frac{1.5-2}{1-0}\right)/1 = -\frac{1}{2}$
\end{itemize}

El sistem queda:

\begin{equation*}
\begin{pmatrix}
	2 & 1 & 0 & 0  \\
	1/2 & 2 & 1/2 & 0  \\
	0 & 1/2 & 2 & 1/2  \\
	0 & 0 & 1 & 1/2
\end{pmatrix}
\begin{pmatrix}
	M_0  \\
	M_1  \\
	M_2  \\
	M_3
\end{pmatrix}
=
\begin{pmatrix}
	6\cdot\frac{3}{10}  \\
	6\cdot\frac{1}{2}  \\
	6\cdot(-1)  \\
	6\cdot(-\frac{1}{2})
\end{pmatrix}
\end{equation*}

Del que obtenemos la solución $M_0=-0.36, M_1=2.52, M_2=-3.72 \text{ y } M_3=0.36$.
Calculamos los trozos finalmente aplicando la fórmula:

$$S_1(x)= M_0\frac{(x_1-x)^3}{6} + M_1\frac{(x-x_0)}{6} + (y_0-\frac{M_0}{6})\frac{x_1-x}{1} + (y_1-\frac{M_1}{6})\frac{x-x_0}{1} = 0.48x^3 - 0.18x^2 + 0.2x$$

Equivalentemente para $C_2$ y $C_3$, obtenemos la solución:

\begin{equation*}
 S(x) =
  \begin{cases}
   0.48x^3 - 0.18x^2 + 0.2x & 							   \text{si } 0 \leq x \leq 1  \\
   -1.04x(x-1)^3 + 1.26(x-1)^2 + 1.28(x-1) - 0.5   &  \text{si } 1 < x \leq 2  \\
   0.68(x-2)^3 - 1.86(x-2)^2 + 0.68(x-2) + 2       &  \text{si } 2 < x \leq 3  \\
  \end{cases}
\end{equation*}
\end{solucion}



\pagebreak

# Implementación en ordenador: Octave

## Spline Lineal

La implementación de la función que nos permite calcular un spline lineal es muy
sencilla:

```octave
function s = SplineLineal(x,y)
  p = diff(y)./diff(x);
  A = [p' y(1:end-1)'];
  s = mkpp(x,A);
end
```

## Splines cuadráticos

Utilizando el **método global**, podemos definir fácilmente una
función que calcule un spline cuadrático de clase 1:

```octave
function s = SplineCuad(x, y, d_k, k)
  # Número de intervalos
  n = length(x) - 1;

  # 1, x, x^2
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
    p = sol(3:-1:1);

    for l = 2:k
      p += sol(l+2).*[1, -2.*x(l), x(l).^2];
    end

    B(k, :) = polyaffine(p,[-x(k) 1]);
  end

  s = mkpp(x,B);
end
```
Otra implementación posible es calcular el spline **a trozos**:

```octave
function z = SplineCuadLocal(x, y, d_k, k)
	s = zeros(length(x)-1, 3);
   d = d_k;

    #Recorremos todos los nodos de n+1 en adelante:

    for i = (k+1):length(x)
		p = (y(i)-y(i-1))/(x(i)-x(i-1));
		q = (p-d)/(x(i)-x(i-1));
		v = [x(i-1) x(i-1)];
		s(i-1,:) = [0 0 y(i-1)]+[0 d -d*x(i-1)]+q*poly(v);
		d = 2*p-d;
	end
    d = d_k;

    #Recorremos todos los nodos desde n hasta el 1:

    for i = 0:(k-2)
		j = k-i;
		p = (y(j)-y(j-1))/(x(j)-x(j-1));
		q = (d-p)/(x(j)-x(j-1));
		v = [x(j-1) x(j)];
		s(j-1,:) = [0 0 y(j-1)]+[0 p -p*x(j-1)]+q*poly(v);
    end

    for i = 1:length(s)
    	s(i,:) = polyaffine(s(i,:), [-x(i), 1]);
    end
    z = mkpp(x, s);
end
```

\pagebreak

## Splines cúbicos

Para el cálculo de splines cúbicos por medio de la segunda derivada nos hemos
valido de la interpolación de splines lineales y de la función `ppint`, que
realiza la integración de un spline.

### Spline sujeto

La función para el cálculo del **spline sujeto** queda:

```octave
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

  m = A\gamma';
  s = ppint(ppint(SplineLineal(x, m')));
end
```

\pagebreak

### Spline natural

Para el **spline natural** calculamos la matriz, resolvemos el sistema e integramos
añadiendo los datos de las derivadas segundas en los extremos:

```octave
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
    A = diag(twoes) + diag(lambda(1:end-1),1) + diag(mu(2:end),-1);
  end

  m = A \ gamma';
  s = ppint(ppint(SplineLineal(x, [0 m' 0])));
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

\begin{definicion}
Dado un intervalo $[a,b]$ y una partición $\{x_i\}_{i = 0...n} \in \mathscr{P}([a,b])$,
se definen, para $1 \leq i \leq n$:

\begin{itemize}
\item $\displaystyle h_i = x_i - x_{i-1}$
\item $\displaystyle p_i = \frac{y_i-y_{i-1}}{h_i}$
\end{itemize}
\end{definicion}


#Bibliografía
- [Quadratic Interpolatory Splines W.J.Kammerer, G. W. Reddien, and R.S. Varga](www.math.kent.edu/~varga/pub/paper_85.pdf)
- [Cubic Spline Interpolation  - Wikiversity](https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation)
- Análisis numérico (Novena Edición) Richard L. Burden y J. Douglas Faires
