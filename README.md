# Oblicua

**Texto e imágenes en perspectiva sobre cualquier foto, directo en el navegador.**

👉 [oblicua.netlify.app](https://oblicua.netlify.app/)

![Oblicua poniendo un título en perspectiva sobre un cuadro](docs/captura.png)

Nació para poner el título de mis covers de piano sobre el atril en los videos. Ahora sirve para cualquier superficie: subes una foto, marcas las cuatro esquinas y el texto (o una imagen) se adapta a su ángulo.

## Qué hace

- Título y subtítulo opcional, cada uno con su propia fuente, y opción de mayúsculas.
- Modo imagen: pon un logo o cualquier PNG sobre la superficie.
- Exporta la foto completa o solo el contenido en **PNG transparente**, listo para un editor de video.
- Lupa de precisión para mover las esquinas, arrastrar y soltar, y pegar desde el portapapeles.
- Todo corre en el navegador: las fotos nunca salen del dispositivo.

## Cómo funciona

- **Homografía** calculada desde cero a partir de las cuatro esquinas; cada píxel de destino se mapea al contenido con la transformación inversa.
- **Remuestreo bilineal con alfa premultiplicado**, para que los bordes del texto queden limpios.
- **Fusión «multiplicar»**, para que el texto parezca impreso sobre la superficie y no pegado encima.
- Sin dependencias ni build: un solo `index.html` con HTML, CSS y JavaScript.

## Uso local

Abre `index.html` en el navegador, o sírvelo con cualquier servidor estático:

```bash
npx serve .
```

## Licencia

MIT © Giovanni Raffa
