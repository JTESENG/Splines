# Metodos

## Compilar el documento

Hay que utilizar `pandoc`:

- `sudo apt-get install pandoc`
- Seguir [estas](http://pandoc.org/installing.html#linux) instrucciones

## Cómo usar Markdown

### Títulos
```md
# Título
## Segundo título
### ...
```

### Fórmulas

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
