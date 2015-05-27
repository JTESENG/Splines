# Metodos

## A Hacer

- [ ] Splines cuadráticos
	- [X] Introducción
	- [X] Descripción del espacio de splines cuadráticos
	- [ ] Interpolación con splines cuadráticos
   	- [ ] A trozos
   	- [X] Global
	- [ ] Ejemplos
- [ ] Splines cúbicos
	- [X] [Construcción de splines clásicos a partir de la segunda derivada](https://en.wikiversity.org/wiki/Cubic_Spline_Interpolation)
   	- [X] Spline sujeto
   	- [ ] Spline natural
   	- [ ] Spline periódico
	- [X] Propiedades de minimización
	- [ ] Ejemplos
- [ ] Implementación en ordenador
	- [ ] Splines cuadráticos
      - [ ] A trozos
      - [X] Global
	- [ ] Splines sujeto
	- [ ] Spline natural
	- [ ] Spline periódico


## Compilar el documento

Hay que utilizar `pandoc`:

- `sudo apt-get install pandoc`
- Seguir [estas](http://pandoc.org/installing.html#linux) instrucciones
- `pandoc Trabajo.md -o Trabajo.pdf`

## Cómo usar Markdown

### Títulos
```md
# Título
## Segundo título
### ...
<!---Comentarios-->
```

### Fórmulas

Tutorial de latex, muy útil para ver como funcionan entornos tales como tablas,
matrices, o derivados del entorno equation:
[The Not So Short Introduction to LaTeX - Tobi Oetiker](https://tobi.oetiker.ch/lshort/lshort.pdf)

Los símbolos se pueden mirar en [DeteXify](http://detexify.kirelabs.org/classify.html).
- Fórmulas dentro del texto: `$ < Código en LaTeX > $`
- Fórmulas con linea propia: `\[ < Código en LaTeX > \]`
- Fórmulas con número para referenciarlas:
     ```md
        \begin{equation}
          < Cosas >
        \end{equation}
      ```
### Formato

 - **Negrita**: `**Negrita**`
 - *Cursiva*: `*Cursiva*`
