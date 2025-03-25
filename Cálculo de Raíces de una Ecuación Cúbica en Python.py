import numpy as np
import math

def calculate_properties(T, P, Tc, Pc, w, R= 0.14305  ):
    # Calcular Tr
    Tr = T / Tc
    
    # Calcular m dependiendo del factor acéntrico w
    if w <= 0.49:
        m = -0.26992 * w**2 + 1.54226 * w + 0.37464
    else:
        m = -0.016666 * w**3 - 0.164423 * w**2 + 1.48503 * w + 0.379642
    
    # Calcular alfa
    alpha = (1 + m * (1 - math.sqrt(Tr)))**2
    
    # Calcular a y b
    a = (0.457235 * R**2 * Tc**2 * alpha) / Pc
    b = (0.077796 * R * Tc) / Pc
    
    # Calcular A y B
    A = (a * P) / (R**2 * T**2)
    B = (b * P) / (R * T)

    # Coeficientes de la ecuación cúbica de Z
    coef = [1, -(1 - B), (A - 2*B - 3*B**2), -(A*B - B**2 - B**3)]
    
    # Resolver la ecuación cúbica
    Z_roots = np.roots(coef)

    # Calcular volumen específico para cada raíz de Z
    volumes = [Z * R * T / P for Z in Z_roots]

    return Tr, m, alpha, a, b, A, B, Z_roots, volumes

# Ejemplo de uso
T = 360        # Temperatura en K
P = 1542.8   # Presión en kPa
Tc = 407.81  # Temperatura crítica en K
Pc = 3629 # Presión crítica en kPa
w = 0.184 # Factor acéntrico

Tr, m, alpha, a, b, A, B, Z_roots, volumes = calculate_properties(T, P, Tc, Pc, w)

print(f"Tr: {Tr}")
print(f"m: {m}")
print(f"Alpha: {alpha}")
print(f"a: {a} kPa·m^6/kmol^2")
print(f"b: {b} m^3/kmol")
print(f"A: {A}")
print(f"B: {B}")
print(f"Raíces de Z: {Z_roots}")
print(f"Volúmenes específicos: {volumes}")