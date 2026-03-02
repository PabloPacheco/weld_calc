# Manual Teórico: Cálculo de Uniones Soldadas según el Método del Área de Garganta

Este documento describe la teoría implementada en el código `weld_calculator.py` para el análisis de tensiones en uniones soldadas de filete. El método sigue el enfoque clásico de *proyección sobre el área de garganta* presentado en textos como Shigley, y se ha extendido para manejar cordones con orientación arbitraria y estados de carga generales (fuerzas y momentos en 3D). Todas las ecuaciones se presentan en unidades del Sistema Internacional (m, N, Pa).

## 1. Introducción

En una soldadura de filete, la sección resistente se idealiza como el **área de garganta**, un plano inclinado a 45° respecto a las caras del cordón. Sobre esta área actúan tensiones normales ($\sigma$) y cortantes ($\tau$) que se calculan proyectando las fuerzas y momentos resultantes. El código presentado automatiza este proceso para un conjunto de cordones, permitiendo:

- Definir cordones mediante sus puntos extremos y su garganta $h$.
- Calcular propiedades geométricas globales (centroide, momentos de inercia, etc.).
- Aplicar múltiples cargas (fuerzas en cualquier dirección) y obtener las tensiones en puntos a lo largo de los cordones.
- Evaluar el factor de seguridad según los criterios de von Mises o Tresca.
- Visualizar los resultados mediante gráficos vectoriales y de barras.

## 2. Geometría del cordón de soldadura

Cada cordón se representa como un segmento rectilíneo en el plano $XY$ (la coordenada $Z$ se reserva para cargas fuera del plano). El cordón tiene una **garganta** $h$ (en metros) que define el espesor efectivo.

### 2.1. Longitud y área de garganta

La longitud del cordón es la distancia entre sus puntos extremos:

$$ l = \|\mathbf{P}_{\text{fin}} - \mathbf{P}_{\text{inicio}}\| \quad [\text{m}] $$

El área de la garganta (superficie de falla) es un rectángulo de dimensiones $l \times (h/\sqrt{2})$:

$$ A = \frac{h}{\sqrt{2}} \, l \quad [\text{m}^2] $$

### 2.2. Centroide del cordón

El centroide del cordón (punto medio del segmento) es:

$$ \mathbf{C} = \frac{\mathbf{P}_{\text{inicio}} + \mathbf{P}_{\text{fin}}}{2} \quad [\text{m}] $$

### 2.3. Momentos de inercia locales

Respecto a su centroide, el área de garganta (rectángulo) tiene dos momentos de inercia principales:

- **Eje longitudinal** (paralelo al cordón):  
  $$ I_{\text{long}} = \frac{1}{12} \left(\frac{h}{\sqrt{2}}\right) l^3 \quad [\text{m}^4] $$
- **Eje transversal** (perpendicular al cordón en el plano de la garganta):  
  $$ I_{\text{trans}} = \frac{1}{12} l \left(\frac{h}{\sqrt{2}}\right)^3 \quad [\text{m}^4] $$

El **momento polar** del cordón (respecto a su centroide) es la suma de ambos:

$$ J_G^{\text{local}} = I_{\text{long}} + I_{\text{trans}} = \frac{1}{\sqrt{2}} \cdot \frac{1}{12} \, h \, l^3 \quad [\text{m}^4] $$

### 2.4. Transformación a ejes globales

Cuando el cordón forma un ángulo $\theta$ con el eje $X$ global, los momentos de inercia respecto a los ejes globales $X, Y$ (con origen en el centroide del cordón) se obtienen mediante rotación:

$$ 
\begin{aligned}
I_{xx}^{\text{loc}} &= I_{\text{long}} \sin^2\theta + I_{\text{trans}} \cos^2\theta \\
I_{yy}^{\text{loc}} &= I_{\text{long}} \cos^2\theta + I_{\text{trans}} \sin^2\theta \\
I_{xy}^{\text{loc}} &= (I_{\text{long}} - I_{\text{trans}}) \sin\theta \cos\theta
\end{aligned}
\quad [\text{m}^4]
$$

donde $\theta = \arctan\left(\frac{\Delta y}{\Delta x}\right)$ con $\Delta x, \Delta y$ las componentes del vector dirección del cordón.

## 3. Propiedades del conjunto de cordones

Para la unión completa (varios cordones), se calculan propiedades globales utilizando el **teorema de Steiner** (transporte de ejes paralelos).

### 3.1. Centroide global

El centroide del conjunto se obtiene como promedio ponderado por área:

$$ 
\mathbf{C}_G = \frac{\sum_{i} A_i \, \mathbf{C}_i}{\sum_{i} A_i} \quad [\text{m}]
$$

### 3.2. Área total

$$ A_{\text{total}} = \sum_i A_i \quad [\text{m}^2] $$

### 3.3. Momento polar equivalente para torsión

El momento polar total (respecto al centroide global) se calcula sumando el momento local de cada cordón más el término de transporte $A_i \, r_i^2$, donde $r_i$ es la distancia del centroide del cordón al centroide global:

$$ J_G = \sum_i \left( J_{G,i}^{\text{local}} + A_i \, r_i^2 \right) \quad [\text{m}^4] $$

### 3.4. Tensor de inercia para flexión

Los momentos y producto de inercia globales (respecto a los ejes $X, Y$ con origen en $\mathbf{C}_G$) se ensamblan de forma análoga:

$$ 
\begin{aligned}
I_{xx} &= \sum_i \left( I_{xx,i}^{\text{loc}} + A_i \, (y_{C_i} - y_G)^2 \right) \\
I_{yy} &= \sum_i \left( I_{yy,i}^{\text{loc}} + A_i \, (x_{C_i} - x_G)^2 \right) \\
I_{xy} &= \sum_i \left( I_{xy,i}^{\text{loc}} + A_i \, (x_{C_i} - x_G)(y_{C_i} - y_G) \right)
\end{aligned}
\quad [\text{m}^4]
$$

Estos valores definen completamente la respuesta de la unión ante momentos flectores.

## 4. Reducción de cargas al centroide

Dado un conjunto de fuerzas $\mathbf{F}_k$ (en Newtons) aplicadas en puntos $\mathbf{P}_k$, se trasladan al centroide global obteniendo una fuerza resultante y un momento resultante:

$$ 
\mathbf{F}_{\text{total}} = \sum_k \mathbf{F}_k \quad [\text{N}]
$$
$$ 
\mathbf{M}_{\text{total}} = \sum_k \left( (\mathbf{P}_k - \mathbf{C}_G) \times \mathbf{F}_k \right) \quad [\text{N·m}]
$$

El vector $\mathbf{M}_{\text{total}}$ tiene componentes $(M_x, M_y, M_z)$.

## 5. Cálculo de tensiones en la garganta

Para un punto de evaluación $\mathbf{P}$ (perteneciente a algún cordón), se determinan las siguientes tensiones.

### 5.1. Tensión cortante primaria (corte directo)

Debida a las componentes de fuerza en el plano $XY$ ( $F_x, F_y$ ):

$$ 
\boldsymbol{\tau}_1 = -\frac{1}{A_{\text{total}}} \begin{bmatrix} F_x \\ F_y \\ 0 \end{bmatrix} \quad [\text{Pa}]
$$

El signo negativo es convencional (la tensión se opone a la fuerza aplicada).

### 5.2. Tensión cortante por torsión (momento $M_z$)

El momento torsor $M_z$ produce una distribución de tensiones cortantes en el plano $XY$:

$$ 
\boldsymbol{\tau}_2 = \frac{-\mathbf{M}_{\text{torsión}} \times (\mathbf{P} - \mathbf{C}_G)}{J_G}, \quad \text{con } \mathbf{M}_{\text{torsión}} = (0,0,M_z) \quad [\text{Pa}]
$$

La tensión cortante total es:

$$ \boldsymbol{\tau} = \boldsymbol{\tau}_1 + \boldsymbol{\tau}_2 $$

### 5.3. Tensión normal (axial + flexión)

La tensión normal en el punto $\mathbf{P}$ tiene tres contribuciones:

- **Axial** (debida a $F_z$):  
  $$ \sigma_{\text{axial}} = \frac{F_z}{A_{\text{total}}} \quad [\text{Pa}] $$

- **Flexión** alrededor de $X$ e $Y$: cuando los ejes no son principales, se utiliza la fórmula general de flexión asimétrica. Sea $\mathbf{r} = \mathbf{P} - \mathbf{C}_G = (x, y, 0)$. La tensión por flexión es:

$$ 
\sigma_{\text{flexión}} = \frac{ (M_x I_{yy} - M_y I_{xy}) \, y + (M_y I_{xx} - M_x I_{xy}) \, x }{ I_{xx} I_{yy} - I_{xy}^2 } \quad [\text{Pa}]
$$

Esta expresión reduce a $\frac{M_x y}{I_{xx}} + \frac{M_y x}{I_{yy}}$ cuando $I_{xy}=0$ (ejes principales).

La tensión normal total es:

$$ \sigma = \sigma_{\text{axial}} + \sigma_{\text{flexión}} $$

### 5.4. Tensión equivalente de von Mises

Para evaluar el riesgo de fluencia en materiales dúctiles, se combinan las tensiones normal y cortante mediante el criterio de von Mises:

$$ 
\sigma_{\text{VM}} = \sqrt{ \sigma^2 + 3 \|\boldsymbol{\tau}\|^2 } \quad [\text{Pa}]
$$

## 6. Factor de seguridad

El factor de seguridad $n$ se define como el cociente entre la tensión de fluencia del material $S_y$ (Pa) y la tensión equivalente:

$$ n = \frac{S_y}{\sigma_{\text{eq}}} $$

Dos criterios están implementados:

- **von Mises**: $\sigma_{\text{eq}} = \sigma_{\text{VM}}$ (recomendado para materiales dúctiles).
- **Tresca** (máxima cortante): $\sigma_{\text{eq}} = 2 \sqrt{(\sigma/2)^2 + \|\boldsymbol{\tau}\|^2 }$.

El código devuelve tanto el factor de seguridad mínimo (punto más crítico) como una lista con los valores en todos los puntos evaluados.

## 7. Visualización de resultados

El módulo incluye varias funciones gráficas:

- **`plot_vectors`**: dibuja los cordones, el centroide y, sobre los puntos de evaluación, vectores (para $\boldsymbol{\tau}$) o círculos coloreados (para $\sigma$ o $\sigma_{\text{VM}}$). Permite escalar y superponer componentes.
- **`plot_sigma_unfolded`**: genera un gráfico de barras con la tensión normal $\sigma$ (en MPa) a lo largo de la unión "desplegada". Cada cordón aparece en una región separada con su color y etiqueta.
- **`plot_n_unfolded`**: similar al anterior, pero muestra el factor de seguridad $n$. Incluye una línea horizontal en $n=1$ para identificar zonas de falla inminente.

## 8. Ejemplo de aplicación

A continuación se muestra un ejemplo básico (para Jupyter Lab) que define tres cordones (horizontal, vertical y oblicuo), aplica dos cargas y visualiza los resultados.

```python
import numpy as np
import matplotlib.pyplot as plt
import weld_calculator as wc

# Crear cordones con nombres
c1 = wc.Weld(name="Horizontal")
c1.start_point = [0, 0, 0]
c1.end_point = [0.05, 0, 0]
c1.h = 0.005
c1.weld_calcs()
c1.color = "black"

c2 = wc.Weld(name="Vertical")
c2.start_point = [0, 0, 0]
c2.end_point = [0, 0.05, 0]
c2.weld_calcs()
c2.color = "green"

c3 = wc.Weld(name="Oblicuo")
c3.start_point = [0.05, 0, 0]
c3.end_point = [0.1, 0.05, 0]
c3.weld_calcs()
c3.color = "blue"

# Unión
wj = wc.Weld_joint()
wj.welds = [c1, c2, c3]
wj.weld_calcs()

# Cargas
F1 = np.array([0, -1000, 0])
p1 = np.array([0.15, 0.05, 0])
F2 = np.array([0, 0, 500])
p2 = np.array([0.025, 0.025, 0])
forces = [F1, F2]
points_F = [p1, p2]

# Puntos de evaluación
pts1 = wc.equidistant_points(c1, 15)
pts2 = wc.equidistant_points(c2, 15)
pts3 = wc.equidistant_points(c3, 15)
weld_points = wc.remove_repeats(np.vstack((pts1, pts2, pts3)))

# Evaluar
wj.weld_evaluation(forces, points_F, weld_points)

# Visualizar tensiones cortantes
wj.plot_vectors(weld_points, plot_type='tau', scale_vectors=0.5, scale_magnitud=1e-6)

# Visualizar tensión normal desplegada
wj.plot_sigma_unfolded(weld_points)

# Visualizar factor de seguridad (suponiendo fluencia de 250 MPa)
wj.plot_n_unfolded(weld_points, fluencia=250e6, criterio='vm')
```

## 9. Referencias

- Shigley, J. E., & Mischke, C. R. (2011). *Diseño en ingeniería mecánica*. McGraw-Hill.
- Blodgett, O. W. (1966). *Design of Welded Structures*. The James F. Lincoln Arc Welding Foundation.


