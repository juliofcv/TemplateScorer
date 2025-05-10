# TemplateScorer

**TemplateScorer** es una clase Java que calcula un **score** de calidad para plantillas de huellas dactilares en formato ISO 19794. Combina técnicas de estadística robusta, análisis multivariante y métricas espaciales para producir un único valor que refleja la consistencia geométrica y direccional de las minucias.

---

## Contenido

- [Descripción](#descripción)  
- [Dependencias](#dependencias)  
- [Uso](#uso)  
- [Detalles de implementación](#detalles-de-implementación)  
  - [Constantes](#constantes)  
  - [scoreTemplate(ReadISO19794 template)](#scoretemplate-readiso19794-template)  
  - [calculateWeights(List<Minutiae> raw)](#calculateweightslistminutiae-raw)  
  - [extractAngles & extractPoints](#extractangles--extractpoints)  
  - [calculateRobustCentroid](#calculaterobustcentroid)  
  - [validateResolution](#validateresolution)  
  - [normalizePoints & calculateScaleFactor](#normalizepoints--calculatescalefactor)  
  - [computePCAAlignment](#computepcaalignment)  
  - [rotatePoints](#rotatepoints)  
  - [geometricMedian](#geometricmedian)  
  - [covariance](#covariance)  
  - [calculateMahalanobisDistances](#calculatemahalanobisdistances)  
  - [calculateTrimmedMean](#calculatetrimmedmean)  
  - [calculateQuadrantScore](#calculatequadrantscore)  
  - [calculateOrientationEntropy](#calculateorientationentropy)  

---

## Descripción

El flujo general de la puntuación es:

1. **Lectura de minucias**  
2. **Pesado** según calidad  
3. **Cálculo de centro robusto** (mediana geométrica)  
4. **Normalización espacial** y **alineación PCA**  
5. **Cálculo de distancias de Mahalanobis**  
6. **Media recortada** de distancias  
7. **Ponderación** por entropía de orientación y distribución espacial (quadrants)

---

## Dependencias

- **Java 8+**  
- `com.gd.bioutils.read.ReadISO19794` (lectura de plantillas ISO 19794)  
- `com.gd.bioutils.entity.Minutiae`  

---

## Uso

```java
import com.gd.bioutils.read.ReadISO19794;
import com.gd.bioutils.scoring.TemplateScorer;

public class Main {
    public static void main(String[] args) {
        ReadISO19794 tpl = ReadISO19794.fromFile("huella.iso19794");
        double score = TemplateScorer.scoreTemplate(tpl);
        System.out.println("Template score: " + score);
    }
}
```

---

## Detalles de implementación

### Constantes

| Constante                    | Descripción                                                                                  |
|------------------------------|----------------------------------------------------------------------------------------------|
| `DEFAULT_RESOLUTION = 500.0` | Resolución por defecto (ppi) si `xRes` o `yRes` están fuera de umbrales.                    |
| `MIN_RESOLUTION_THRESHOLD`   | Mínimo valor válido de resolución (100 ppi).                                                 |
| `MAX_RESOLUTION_THRESHOLD`   | Máximo valor válido de resolución (10000 ppi).                                               |
| `MIN_QUALITY_WEIGHT = 0.01`  | Peso mínimo para una minucia (si `quality/100` es menor).                                    |
| `TRIM_RATIO = 0.1`           | Fracción de datos a recortar en el cálculo de la media (10 %).                               |
| `HIST_BINS = 12`             | Número de contenedores para el histograma de ángulos.                                        |
| `BETA = 0.5`                 | Factor de ponderación para la entropía de orientación.                                       |
| `MAX_MEDIAN_ITER = 10`       | Iteraciones máximas para el algoritmo de mediana geométrica.                                 |
| `MEDIAN_TOL = 1e-6`          | Tolerancia de convergencia para la mediana geométrica.                                       |

### `scoreTemplate(ReadISO19794 template)`

Procesa una plantilla completa:

1. Extrae listas de **minutias**, **pesos**, **ángulos** y **coordenadas**.  
2. Calcula un **centro robusto** (mediana geométrica) con `geometricMedian`[^1].  
3. Valida y normaliza la resolución (`validateResolution`).  
4. Centra y **normaliza** puntos en escala (`normalizePoints`).  
5. Alinea con **PCA** (`computePCAAlignment`)[^2][^3].  
6. Rota los puntos según el ángulo PCA (`rotatePoints`).  
7. Calcula distancias de **Mahalanobis** (`calculateMahalanobisDistances`)[^4].  
8. Obtiene la **media recortada** de las distancias (`calculateTrimmedMean`)[^5].  
9. Calcula **quadrantScore** (distribución espacial).  
10. Calcula **orientationEntropy** (entropía de ángulos) (`calculateOrientationEntropy`)[^6].  
11. Combina como  
    ```
    score = globalScore * (1 + BETA * orientationEntropy) * quadrantScore
    ```  

### `calculateWeights(List<Minutiae> raw)`

Convierte la calidad de cada minutia (0–100) en un peso en [MIN_QUALITY_WEIGHT, 1].

### `extractAngles` & `extractPoints`

Lee ángulo y coordenadas `(x,y)` de cada minutia.

### `calculateRobustCentroid`

Invoca `geometricMedian(points, weights)` para hallar el punto que **minimiza la suma ponderada de distancias**[^1].

### `validateResolution`

Garantiza que `xRes` y `yRes` estén entre `MIN_RESOLUTION_THRESHOLD` y `MAX_RESOLUTION_THRESHOLD`.

### `normalizePoints` & `calculateScaleFactor`

- Calcula un **factor de escala** robusto en base a distancias al centro y pesos.  
- Traslada y divide por resolución y escala para estandarizar.

### `computePCAAlignment`

Calcula el **primer componente principal** mediante la matriz de covarianza[^2][^3]:
```java
covXX = Σ x² / n; covXY = Σ x·y / n; covYY = Σ y² / n;
angle = 0.5 * atan2(2·covXY, covXX - covYY);
```

### `rotatePoints`

Rota cada punto `(x,y)` por `angle`:
```
x' = x·cos – y·sin;  y' = x·sin + y·cos
```

### `geometricMedian`

Implementa el **algoritmo de Weiszfeld** (mediana geométrica ponderada)[^1]:
```java
for iter hasta MAX_MEDIAN_ITER:
  inv = (d < MEDIAN_TOL) ? 1 : weight / d;
  next = Σ (p·inv) / Σ inv;
  si ‖next–median‖ < MEDIAN_TOL ⇒ break;
```

### `covariance`

Matriz de covarianza ponderada:
```
c00 = Σ w·dx² / Σw;  c01 = Σ w·dx·dy / Σw;  c11 = Σ w·dy² / Σw;
```

### `calculateMahalanobisDistances`

- Invierte la matriz de covarianza y calcula  
  `d² = [x y]·Cov⁻¹·[x; y]`  
- Distancia = √max(d², 0)  
- Ordena ascendentemente.

### `calculateTrimmedMean`

- Ordena distancias  
- Elimina `⌊n·TRIM_RATIO⌋` de extremos  
- Devuelve media de las restantes (estadística robusta)[^5].

### `calculateQuadrantScore`

Divide el bounding box en `k×k` subcuadrantes (`k=3`), calcula la media recortada de distancias al origen en cada bloque, y devuelve la raíz geométrica de sus productos.

### `calculateOrientationEntropy`

- Cuenta frecuencias de ángulos en `HIST_BINS` contenedores  
- Calcula entropía de Shannon[^6]:  
  `–Σ₁ᵇ pᵢ·log₂(pᵢ)`  
- Normaliza por `log₂(HIST_BINS)` para obtener [0,1].

---

## Referencias

[^1]: Weiszfeld, E. (1937). *Sur le point pour lequel la somme des distances de n points donnés est minimum.*  
[^2]: Pearson, K. (1901). *On lines and planes of closest fit to systems of points in space.*  
[^3]: Hotelling, H. (1933). *Analysis of a complex of statistical variables into principal components.*  
[^4]: Mahalanobis, P.C. (1936). *On the generalized distance in statistics.*  
[^5]: Huber, P.J. (1981). *Robust Statistics.*  
[^6]: Shannon, C.E. (1948). *A Mathematical Theory of Communication.*  
