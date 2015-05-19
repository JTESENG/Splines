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

# Implementación en ordenador
