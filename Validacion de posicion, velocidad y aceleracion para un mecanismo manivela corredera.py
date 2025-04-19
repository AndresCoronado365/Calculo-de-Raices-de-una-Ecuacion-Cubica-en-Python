import numpy as np

def newton_raphson_4bar(L1y, L2, L3, theta2_deg, theta_semilla, tol=1e-10, max_iter=100):

    theta2 = np.radians(theta2_deg)
    theta3 = np.radians(theta_semilla[0])
    L1x = theta_semilla[1]
    
    for _ in range(max_iter):
        F1 = L2 * np.cos(theta2) + L3 * np.cos(theta3) - L1x
        F2 = L2 * np.sin(theta2) + L3 * np.sin(theta3) + L1y
        
        J = np.array([
            [-L3 * np.sin(theta3), -1],
            [L3 * np.cos(theta3), 0]
        ])
        
        F = np.array([F1, F2])
        delta = np.linalg.solve(J, F)
        
        theta3 -= delta[0]
        L1x -= delta[1]
        
        if np.linalg.norm(delta) < tol:
            break
    
    return np.degrees(theta3), L1x

def calcular_velocidades(L1y, L2, L3, theta2_deg, theta3_deg, L1x, theta2_dot):
   
    theta2 = np.radians(theta2_deg)
    theta3 = np.radians(theta3_deg)
    
    A = np.array([
        [-L3 * np.sin(theta3), -1],
        [L3 * np.cos(theta3), 0]
    ])
    
    b = np.array([
        L2 * np.sin(theta2) * theta2_dot,
        -L2 * np.cos(theta2) * theta2_dot
    ])
    
    # Resolver el sistema A * [θ̇3, L1ẋ] = b
    derivatives = np.linalg.solve(A, b)
    
    return derivatives[0], derivatives[1]  # θ̇3, L1ẋ


def calcular_aceleraciones(L1y, L2, L3, theta2_deg, theta3_deg, L1x, 
                          theta2_dot, theta3_dot, L1x_dot, theta2_ddot=0):
    """
    Calcula las aceleraciones angulares (θ̈3) y aceleración lineal (L1ẍ) del mecanismo.
    
    Retorna:
    - theta3_ddot en rad/s², L1x_ddot en mm/s²
    """
    theta2 = np.radians(theta2_deg)
    theta3 = np.radians(theta3_deg)
    
    # Matriz de coeficientes A (la misma que para velocidades)
    A = np.array([
        [-L3 * np.sin(theta3), -1],
        [L3 * np.cos(theta3), 0]
    ])
    
    # Términos conocidos del lado derecho
    B = np.array([
        L3 * np.cos(theta3) * theta3_dot**2,
        L3 * np.sin(theta3) * theta3_dot**2
    ])
    
    C = np.array([
        L2 * np.sin(theta2) * theta2_ddot + L2 * np.cos(theta2) * theta2_dot**2,
        -L2 * np.cos(theta2) * theta2_ddot + L2 * np.sin(theta2) * theta2_dot**2
    ])
    
    # Resolver el sistema A * [θ̈3, L1ẍ] = B + C
    second_derivatives = np.linalg.solve(A, B + C)
    
    return second_derivatives[0], second_derivatives[1]  # θ̈3, L1ẍ


# Parámetros del mecanismo
L1y = 0 # Longitud del eslabón fijo en y (mm)
L2 = 90 # Eslabón de entrada (mm)
L3 = 106 # Eslabón de acoplamiento (mm)

theta2_dot = 2 * np.pi *20/ 60  # Velocidad angular de entrada en rad/s (20 RPM)
theta2_ddot = 0  # Se asume constante

# Ángulos de entrada y valores semilla
theta2_deg_list = [-35.19512195



 ,-168



,-59]
theta_semilla = [20,50]  # [theta3 inicial en grados, L1x inicial en mm]

# Cálculo para cada posición
for theta2_deg in theta2_deg_list:
    theta3_sol, L1x_sol = newton_raphson_4bar(L1y, L2, L3, theta2_deg, theta_semilla)
    
    theta3_dot, L1x_dot = calcular_velocidades(L1y, L2, L3, theta2_deg, theta3_sol, L1x_sol, theta2_dot)
    theta3_ddot, L1x_ddot = calcular_aceleraciones(L1y, L2, L3, theta2_deg, theta3_sol, L1x_sol, 
                                                  theta2_dot, theta3_dot, L1x_dot, theta2_ddot)
    
    print(f"\nPara θ2 = {theta2_deg}°:")
    print(f"  θ3 = {theta3_sol:.4f}°, L1x = {L1x_sol:.4f} mm")
    print(f"  θ̇3 = {theta3_dot:.4f} rad/s, L1ẋ = {L1x_dot:.4f} mm/s")
    print(f"  θ̈3 = {theta3_ddot:.4f} rad/s², L1ẍ = {L1x_ddot:.4f} mm/s²")
