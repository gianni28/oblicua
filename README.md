# Oblicua

**Texto e imágenes en perspectiva sobre cualquier foto o video, directo en el navegador.**

👉 [oblicua.netlify.app](https://oblicua.netlify.app/)

![Oblicua poniendo un título en perspectiva sobre un cuadro](docs/captura.png)

Nació para poner el título de mis covers de piano sobre el atril en los videos. Ahora sirve para cualquier superficie: subes una foto o un video, marcas las cuatro esquinas y el texto (o una imagen) se adapta a su ángulo.

## Qué hace

- Título y subtítulo opcional, cada uno con su propia fuente, y opción de mayúsculas.
- Modo imagen: pon un logo o cualquier PNG sobre la superficie.
- **Videos**: eliges un fotograma, ubicas el letrero y descargas el video completo en MP4 con el letrero puesto y el audio original. Puedes decidir en qué segundo aparece y desaparece, con transición suave.
- Exporta la foto completa o solo el contenido en **PNG transparente**, listo para un editor de video.
- Lupa de precisión para mover las esquinas, arrastrar y soltar, y pegar desde el portapapeles.
- Todo corre en el navegador: las fotos y videos nunca salen del dispositivo.

## Cómo funciona

- **Homografía** calculada desde cero a partir de las cuatro esquinas; cada píxel de destino se mapea al contenido con la transformación inversa.
- **Remuestreo bilineal con alfa premultiplicado**, para que los bordes del texto queden limpios.
- **Fusión «multiplicar»**, para que el texto parezca impreso sobre la superficie y no pegado encima.
- **Video con WebCodecs**: el contenido deformado se calcula una sola vez como capa transparente y se compone sobre cada fotograma en canvas; [Mediabunny](https://mediabunny.dev) se encarga de leer el archivo, decodificar, codificar a MP4 (H.264 cuando el navegador lo permite) y copiar el audio sin recomprimirlo. Respeta la rotación de los videos grabados en vertical.
- Sin build: un `index.html` con HTML, CSS y JavaScript, más Mediabunny en `vendor/`, que solo se descarga al exportar un video.

## Uso local

Sírvelo con cualquier servidor estático (abrir el archivo directamente no permite cargar el módulo de video):

```bash
npx serve .
```

## Licencia

MIT © Giovanni Raffa. Mediabunny se distribuye bajo MPL-2.0 (`vendor/MEDIABUNNY-LICENSE`).
