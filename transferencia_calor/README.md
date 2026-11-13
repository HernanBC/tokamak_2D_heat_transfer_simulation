# Transferencia de calor, materiales y refrigeración
## 1. Descripción

&emsp;En esta sección van a estar las propiedades termofísicas, las ecuaciones de transferencia de calor y condiciones de refrigeración que necesitemos para desarrollar el modelo térmico 2D del tokamak con parámetros de referencia del ITER.

&emsp;La información se organiza para que luego se implemente en el código, con datos de ITER Organization, artículos científicos que encontré, el NRL Plasma Formulary y bibliografía de heat transferencie.

&emsp;Para la primera simulación como vamos a considerar un modelo estacionario, sin movimiento macroscópico del plasma, con calentamiento óhmico y pérdidas radiativas por bremsstrahlung. Se va a emplear una temperatura común para los electrones e iones (quedando pendiente definir si varía espacialmente)

&emsp;Las ecuaciones incluidas tienen que aplicarse dentro de las condiciones válidas que vienen indicadas en sus respectivas fuentes.

## 2. Conductividad térmica del plasma

### 2.1. Transporte térmico anisotrópico

&emsp;En un plasma magnetizado, el transporte térmico depende de la orientación del gradiente de temperatura con respecto al campo magnético.

Entonces se distinguen dos conductividades principales:

- k∥: conductividad paralela al campo magnético [W/(m·K)].
- k⊥: conductividad perpendicular al campo magnético [W/(m·K)].

&emsp;La teoría clásica de Braginskii lo que hace es decribir el transporte colisional mediante coeficientes que dependen de las propiedades de los electrones e iones y de la magnetización del plasma.

&emsp;Las expresiones generales que hay documentadas en PlasmaPy son:

**Conductividad electrónica:**

kₑ = k̂ₑ · (nₑ · kB² · Tₑ · τₑ) / mₑ

**Conductividad iónica:**

kᵢ = k̂ᵢ · (nᵢ · kB² · Tᵢ · τᵢ) / mᵢ

Donde:


| Variable | Descripción | Unidad |
|---|---|---|
| nₑ | Densidad numérica de electrones | m⁻³ |
| nᵢ | Densidad numérica de iones | m⁻³ |
| Tₑ | Temperatura electrónica | K |
| Tᵢ | Temperatura iónica | K |
| kB | Constante de Boltzmann | J/K |
| mₑ | Masa del electrón | kg |
| mᵢ | Masa del ion | kg |
| τₑ | Tiempo de colisión electrónica | s |
| τᵢ | Tiempo de colisión iónica | s |
| k̂ | Coeficiente adimensional de transporte | — |


Y también es importante que los coeficientes adimensionales se selecciondn según la especie, el modelo colisional y la dirección del campo magnético.

&emsp;La ley de conducción térmica anisotrópica se puede expresar mediante:

**q'' = -K · grad(T)**

Ahora, para el transporte conductivo paralelo y perpendicular al campo:

**K = k_perp · I + (k_par - k_perp) · (b ⊗ b)**

donde - q'': vector de flujo térmico [W/m²].
- K: tensor de conductividad térmica [W/(m·K)].
- I: tensor identidad.
- b: vector unitario en la dirección del campo magnético.
- ⊗: producto tensorial.

&emsp;Acá se considera la parte simétrica del transporte conductivo y no incluye en explícito el término transversal de Braginskii.

**Condiciones y limitaciones:**

- El modelo supone un plasma completamente ionizado y transporte colisional.
- No representa por sí solo el transporte turbulento característico de un tokamak.
- Si se considera una temperatura estrictamente uniforme en todo el plasma, entonces el gradiente térmico sería nulo y la conducción no podría evaluarse como mecanismo de transporte espacial.
- Los coeficientes numéricos específicos de k∥ y k⊥ hacen falta de seleccionarse y comprobarse para las condiciones del proyecto.

**Referencias**

PlasmaPy. *Classical transport theory (Braginskii)*.

https://docs.plasmapy.org/en/stable/formulary/braginskii.html

PlasmaPy. *ClassicalTransport — API Reference*.

https://docs.plasmapy.org/en/stable/api/plasmapy.formulary.braginskii.ClassicalTransport.html

## 3. Materiales de la primera pared

### 3.1. Configuración de materiales

De acuerdo con ITER Organization y Sardo et al. (2026), la primera pared actualizada utiliza los siguientes materiales:

| Material | Función |
|---|---|
| Tungsteno (W) | Superficie expuesta al plasma |
| CuCrZr-IG | Disipador térmico con canales de refrigeración |
| Acero inoxidable 316L(N)-IG | Estructura de soporte |

Luego en 2023, ITER decidió sustituir el berilio por tungsteno como material de revestimiento de la primera pared.

### 3.2. Conductividad térmica del tungsteno

Los siguientes valores fueron tomados de la tabla 1 de Wang et al. (2020).

| Temperatura (°C) | Conductividad [W/(m·K)] |
|---:|---:|
| 20 | 173 |
| 200 | 156 |
| 500 | 133 |
| 800 | 118 |
| 1000 | 111 |
| 1500 | 101 |
| 1800 | 99 |

Intervalo tabulado: 20–1800 °C.

### 3.3. Conductividad térmica del CuCrZr

Datos de Wang et al. (2020), tabla 1.

| Temperatura (°C) | Conductividad [W/(m·K)] |
|---:|---:|
| 20 | 326 |
| 250 | 343 |
| 500 | 348 |

Intervalo tabulado: 20–500 °C.

Estas propiedades de ambas tablas son de los materiales utilizados en el estudio térmico citado.
### 3.4. Acero inoxidable 316L(N)-IG

Se puede utilizar en su estructura en la configuración actualizada de los paneles de primera pared.

La tabla específica de conductividad térmica en función de la temperatura está pendiente de validar para esta variante del material.

### 3.5. Conducción térmica en las paredes

Para sólidos isotrópicos se utiliza la ley de Fourier:

**q'' = -k · grad(T)**

Para conducción unidimensional:

**Q̇ = -k · A · (dT/dx)**

Cuando la conductividad depende de la temperatura, puede utilizarse su promedio integral en el intervalo, siempre que se consider:

**k_prom = [1 / (T₂ - T₁)] · ∫₍T₁₎⁽T₂⁾ k(T) dT**

**Referencias:**

ITER Organization. *Blanket*.

https://www.iter.org/machine/blanket

Sardo, M., et al. (2026). *EU ITER first wall panel CuCrZr material procurement and characterization*. Fusion Engineering and Design, 222, 115417.

https://doi.org/10.1016/j.fusengdes.2025.115417

Wang, S., et al. (2020). *Thermal damage of tungsten-armored plasma-facing components under high heat flux loads*. Scientific Reports, 10, 1359. Tabla 1.

https://doi.org/10.1038/s41598-020-57852-8

Çengel, Y. A., & Ghajar, A. J. (2020). *Heat and Mass Transfer: Fundamentals & Applications* (6th ed.). McGraw-Hill. Capítulos 2 y 3.

## 4. Pérdidas radiativas del plasma

### 4.1. Radiación bremsstrahlung

La radiación de frenado ocurre por la interacción de los electrones con los campos eléctricos de los iones.

El NRL Plasma Formulary presenta la siguiente ecuación para la potencia radiada por unidad de volumen:


**P_Br = 1.69 × 10⁻³² · Nₑ · √Tₑ · Σ[Z² · N(Z)]**

donde las densidades numéricas van en cm⁻³, l temperatura electrónica en eV y la potencia radiada por unidad de volumen en W/cm³.

Ahora, expresada con densidades en m⁻³ y potencia en W/m³:

La potencia radiada por unidad de volumen, expresada en unidades del SI:

**S_rad = 1.69 × 10⁻³⁸ · nₑ · √(Tₑ,eV) · Σ[Z² · n_Z]**

Bien, para el caso particular de deuterio completamente ionizado (Z = 1), sin impurezas y con nᵢ = nₑ = n:

**S_rad = 1.69 × 10⁻³⁸ · n² · √(T_eV)**
La temperatura tiene que ir en eV para utilizar esta constante y la potencia radiada total se obtiene integrando la fuente volumétrica:

**P_rad = ∫(Vₚ) S_rad · dV**

Nota: en el balance energético del plasma, la radiación se considera una pérdida.

Parte de la radiación que emite el plasma llega a las paredes y a otrso componentes de la cámara, entonces la cantidad de calor que recibe la superficie depende de su ubicación, geometría y de cómo se distribuye la radiación dentro del plasma.
**Como limitaciones:** esta aproximación no incluye radiación de líneas de impurezas, recombinación, radiación ciclotrónica ni efectos de absorción dentro del plasma.

**Referencia:**

Naval Research Laboratory. (2023). *NRL Plasma Formulary*. Sección *Radiation*, ecuación (30).

https://www.nrl.navy.mil/News-Media/Publications/NRL-PlasmaFormulary/

## 5. Transferencia de calor hacia la primera pared y el divertor

El calor que sale del plasma hacia las paredes y otros componenetes no se puede calcular como si atravesara un material sólid, porque el plasma está confinado por un campo magnético que influye en la forma en que se trasnfiere la energía.
El programa SOLPS-ITER permite estudiar cmo se trasnforma la energía en la zona externa del plasma, y para esto se considera el flujo de enrgía de los electrones e iones hacia las superficies del divertor, incluyendo lo que ocurre en la región cercana a estas superficies.
Para la simulación simplificada, el objetivo es obtener las cargas térmicas superficiales:

- **q''_pared:** flujo de calor que llega a la primera pared [W/m²].
- **q''_divertor:** flujo de calor que llega al divertor [W/m²].
Una vez que ya se conoce el flujo de calor que llega a una superficie, este se puede utilizar para calcular la temperatura de la pared.
Bien, igual hace falta definir cuanto calor va a recibir la primera pared y cuanto llega al divertor. Para poder establecerla, necesitamos encontrar alguna referencia o así que lo respalde o que acordemos algún valor

También hay que evitar contabilizar dos veces la energía radiada cuando se incluyan simultáneamente las pérdidas volumétricas del plasma y la potencia absorbida por las paredes.

**Referencias:**

SOLPS-ITER Tutorials. *Adjusting SOLPS-ITER input*.

https://solps-tutorials.readthedocs.io/my_first_simulation/Adjusting_SOLPS-ITER_input/

SOLPS-ITER Tutorials. *Boundary-condition parameters*.

https://solps-tutorials.readthedocs.io/extras/b2input/develop/b2.parameters.html

## 6. Refrigeración de las paredes

### 6.1. Condiciones de entrada del refrigerante

ITER Organization indica que el sistema de refrigeración del blanket utiliza agua con las siguientes condiciones de entrada:

| Parámetro | Valor | Unidad SI |
|---|---:|---|
| Temperatura de entrada | 70 °C | 343.15 K |
| Presión de entrada | 4 MPa | 4.0 × 10⁶ Pa |

Estos valores corresponden a los del agua utilizada para enfrair el blanket del ITER, entonces no necesariamente son los mismos que se utlizan en el sistema de refigeración del divertor.

### 6.2. Propiedades necesarias

Para determinar la transferencia de calor hacia el refrigerante se necesitan:

| Propiedad | Símbolo | Unidad |
|---|---|---|
| Densidad | ρ_f | kg/m³ |
| Viscosidad dinámica | μ_f | Pa·s |
| Conductividad térmica | k_f | W/(m·K) |
| Calor específico | cp_f | J/(kg·K) |
| Caudal másico | ṁ | kg/s |

Las propiades del agua deben obtenerse considerando la temperatura y presión a las que trabaja el sistema de refrigeración.
Por ahora, falta consultar en la base de datos NIST los valores de densidad, viscosidad dinámica, conudctividad térmica y calor específico del agua a 4MPa, así como dterminar el caudal másico que circula por los canales que se van a modelar.

### 6.3. Convección forzada interna

Las relaciones utilizadas para determinar el coeficiente de transferencia de calor por convección son:

**Número de Reynolds:**

Re = (ρ_f · v · D_h) / μ_f

**Número de Prandtl:**

Pr = (μ_f · cp_f) / k_f

**Número de Nusselt:**

Nu = (h · D_h) / k_f

Despejando el coeficiente de convección:

**h = (Nu · k_f) / D_h**

Donde:

- **ρ_f:** densidad del refrigerante [kg/m³].
- **v:** velocidad media del refrigerante [m/s].
- **D_h:** diámetro hidráulico del canal [m].
- **μ_f:** viscosidad dinámica del refrigerante [Pa·s].
- **cp_f:** calor específico del refrigerante [J/(kg·K)].
- **k_f:** conductividad térmica del refrigerante [W/(m·K)].
- **h:** coeficiente de convección [W/(m²·K)].

La correlación para calcular Nu se va a seleccionar de acuerdo al régimen del flujo, la geometría de los canales y las condiciones de trasnferencia de calor.
Para el balance de energía del agua se considera un fluko estacionario, despreciando  los cambios de enrgía cinética y potencial, y el trabajo.

Q̇_f = ṁ · (h_salida - h_entrada)

Para un calor específico aproximadamente constante:

**Q̇_f = ṁ · cp_f · (T_salida - T_entrada)**
**Referencias:**

ITER Organization. *Blanket*.

https://www.iter.org/machine/blanket

ITER Organization. *Cooling Water*.

https://www.iter.org/machine/supporting-systems/cooling-water

National Institute of Standards and Technology. *Thermophysical Properties of Fluid Systems*.

https://webbook.nist.gov/chemistry/fluid/

Çengel, Y. A., & Ghajar, A. J. (2020). *Heat and Mass Transfer: Fundamentals & Applications* (6th ed.). McGraw-Hill. Capítulos 6 y 8.


