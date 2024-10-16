import numpy as np
import matplotlib.pyplot as plt

# Constantes et conditions initiales
m1, m2, l1, l2 = 1, 1, 1, 1
Y0 = [np.pi/4, np.pi/4, 0, 0]  # Conditions initiales [theta1, theta2, dtheta1, dtheta2]
T = 2  # Durée totale de la simulation
dt = 0.01  # Pas de temps
N = int(T / dt)  # Nombre de pas

# Fonction définissant le système de pendules
def F(Y, l1, l2, m1, m2):
    th1, th2, dth1, dth2 = Y
    g = 9.81  # Gravité
    d2th1 = (-dth1**2*l1*m2*np.sin(2*th1 - 2*th2)/2 - dth2**2*l2*m2*np.sin(th1 - th2) - g*m1*np.sin(th1) - g*m2*np.sin(th1)/2 - g*m2*np.sin(th1 - 2*th2)/2)/(l1*(m1 - m2*np.cos(th1 - th2)**2 + m2))
    d2th2 = (dth1**2*l1*m1*np.sin(th1 - th2) + dth1**2*l1*m2*np.sin(th1 - th2) + dth2**2*l2*m2*np.sin(2*th1 - 2*th2)/2 - g*m1*np.sin(th2)/2 + g*m1*np.sin(2*th1 - th2)/2 - g*m2*np.sin(th2)/2 + g*m2*np.sin(2*th1 - th2)/2)/(l2*(m1 - m2*np.cos(th1 - th2)**2 + m2))
    return np.array([dth1, dth2, d2th1, d2th2])

# Méthodes de résolution numérique
def methode_euler(F, Y0, dt, N):
    Y = np.zeros((4, N))
    Y[:, 0] = Y0
    for i in range(1, N):
        Y[:, i] = Y[:, i - 1] + dt * F(Y[:, i - 1], l1, l2, m1, m2)
    return Y

def runge_kutta_4(F, Y0, dt, N):
    Y = np.zeros((4, N))
    Y[:, 0] = Y0
    for i in range(1, N):
        k1 = F(Y[:, i - 1], l1, l2, m1, m2)
        k2 = F(Y[:, i - 1] + 0.5 * dt * k1, l1, l2, m1, m2)
        k3 = F(Y[:, i - 1] + 0.5 * dt * k2, l1, l2, m1, m2)
        k4 = F(Y[:, i - 1] + dt * k3, l1, l2, m1, m2)
        Y[:, i] = Y[:, i - 1] + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)
    return Y

def verlet(F, Y0, dt, N):
    Y = np.zeros((4, N))
    Y[:, 0] = Y0
    accelerations = F(Y[:, 0], l1, l2, m1, m2)[2:]
    Y[0:2, 1] = Y[0:2, 0] + dt * Y[2:, 0]
    Y[2:, 1] = Y[2:, 0] + dt * accelerations

    for i in range(2, N):
        nouvelles_accelerations = F(Y[:, i - 1], l1, l2, m1, m2)[2:]
        Y[0:2, i] = 2 * Y[0:2, i - 1] - Y[0:2, i - 2] + dt**2 * nouvelles_accelerations
        Y[2:, i] = Y[2:, i - 1] + 0.5 * dt * (nouvelles_accelerations + accelerations)
        accelerations = nouvelles_accelerations 

    return Y

# Calcul des solutions
Y_euler = methode_euler(F, Y0, dt, N)
Y_rk4 = runge_kutta_4(F, Y0, dt, N)
Y_verlet = verlet(F, Y0, dt, N)

# Conversion en coordonnées cartésiennes pour le tracé
x1_euler, y1_euler = l1 * np.sin(Y_euler[0, :]), -l1 * np.cos(Y_euler[0, :])
x2_euler, y2_euler = x1_euler + l2 * np.sin(Y_euler[1, :]), y1_euler - l2 * np.cos(Y_euler[1, :])

x1_rk4, y1_rk4 = l1 * np.sin(Y_rk4[0, :]), -l1 * np.cos(Y_rk4[0, :])
x2_rk4, y2_rk4 = x1_rk4 + l2 * np.sin(Y_rk4[1, :]), y1_rk4 - l2 * np.cos(Y_rk4[1, :])

x1_verlet, y1_verlet = l1 * np.sin(Y_verlet[0, :]), -l1 * np.cos(Y_verlet[0, :])
x2_verlet, y2_verlet = x1_verlet + l2 * np.sin(Y_verlet[1, :]), y1_verlet - l2 * np.cos(Y_verlet[1, :])

# Tracé des trajectoires
plt.figure(figsize=(12, 9))
plt.plot(x1_euler, y1_euler, 'r--', label='Première masse Euler')
plt.plot(x2_euler, y2_euler, 'r-', label='Deuxième masse Euler')
plt.plot(x1_rk4, y1_rk4, 'g--', label='Première masse RK4')
plt.plot(x2_rk4, y2_rk4, 'g-', label='Deuxième masse RK4')
plt.plot(x1_verlet, y1_verlet, 'b--', label='Première masse Verlet')
plt.plot(x2_verlet, y2_verlet, 'b-', label='Deuxième masse Verlet')
plt.xlabel('Position X')
plt.ylabel('Position Y')
plt.title('Trajectoires des deux masses du double pendule (0 à 2s)')
plt.legend()
plt.axis('equal')
plt.savefig('trajectoires_3methodes.png')
plt.show()
