import matplotlib
import tkinter as tk

matplotlib.use('TkAgg')  # Utilisez un backend interactif
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rc('figure', figsize=(12, 9))
plt.rcParams['animation.writer'] = 'ffmpeg'

# Paramètres du système et simulation
l1, l2, m1, m2 = 1, 1, 1, 1
T = 20
dt = T / 1000
alpha_values = [10**(-9), 10**(-7), 0.005, 0.05, 1, 1.2, 1.5, 1.7, 1.9, 2]


def F(Y, l1, l2, m1, m2):
    th1, th2, dth1, dth2 = Y
    g = 9.81  # Constante gravitationnelle

    d2th1 = (-dth1**2*l1*m2*np.sin(2*th1 - 2*th2)/2 - dth2**2*l2*m2*np.sin(th1 - th2) - g*m1*np.sin(th1) - g*m2*np.sin(th1)/2 - g*m2*np.sin(th1 - 2*th2)/2)/(l1*(m1 - m2*np.cos(th1 - th2)**2 + m2))
    d2th2 = (dth1**2*l1*m1*np.sin(th1 - th2) + dth1**2*l1*m2*np.sin(th1 - th2) + dth2**2*l2*m2*np.sin(2*th1 - 2*th2)/2 - g*m1*np.sin(th2)/2 + g*m1*np.sin(2*th1 - th2)/2 - g*m2*np.sin(th2)/2 + g*m2*np.sin(2*th1 - th2)/2)/(l2*(m1 - m2*np.cos(th1 - th2)**2 + m2))
    return np.array([dth1, dth2, d2th1, d2th2])


# Méthode Runge-Kutta d'ordre 4
def method_RK4(Y, l1, l2, m1, m2, dt, T):
    N = int(T / dt)
    temps = np.linspace(0, T, N)
    sol = np.zeros((4, N))
    sol[:, 0] = Y
    for i in range(N - 1):
        k1 = F(sol[:, i], l1, l2, m1, m2)
        k2 = F(sol[:, i] + dt / 2 * k1, l1, l2, m1, m2)
        k3 = F(sol[:, i] + dt / 2 * k2, l1, l2, m1, m2)
        k4 = F(sol[:, i] + dt * k3, l1, l2, m1, m2)
        sol[:, i + 1] = sol[:, i] + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return temps, sol


# Fonction pour créer une animation combinée sur un seul plot
def create_combined_animation(alpha_values, l1, l2, m1, m2, dt, T, filename):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(-(l1 + l2) * 1.1, (l1 + l2) * 1.1)
    ax.set_ylim(-(l1 + l2) * 1.1, (l1 + l2) * 1.1)

    lines = []
    data = []

    # Création de lignes distinctes pour chaque alpha dans un même axe
    colors = plt.cm.viridis(np.linspace(0, 1, len(alpha_values)))  # Utiliser une colormap pour différentes couleurs
    for alpha, color in zip(alpha_values, colors):
        Y0 = [alpha, alpha, 0, 0]
        temps, sol = method_RK4(Y0, l1, l2, m1, m2, dt, T)
        x1 = l1 * np.sin(sol[0, :])
        y1 = -l1 * np.cos(sol[0, :])
        x2 = x1 + l2 * np.sin(sol[1, :])
        y2 = y1 - l2 * np.cos(sol[1, :])
        line, = ax.plot([], [], 'o-', lw=2, color=color, label=f'Alpha = {alpha}')
        lines.append(line)
        data.append((x1, y1, x2, y2))

    ax.legend()

    def init():
        for line in lines:
            line.set_data([], [])
        return lines

    def animate(i):
        for line, (x1, y1, x2, y2) in zip(lines, data):
            line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
            print(f"Enregistrement de l'image {i + 1}/{len(temps)}")
        return lines

    ani = animation.FuncAnimation(fig, animate, init_func=init, frames=len(temps), blit=True)
    ani.save(filename, writer='ffmpeg', fps=30)
    plt.close()

# Fonction pour tracer les trajectoires des pendules
def plot_trajectories(temps, sol, filename, alpha):
    theta1, theta2 = sol[0, :], sol[1, :]
    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)
    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)
    plt.figure()
    plt.plot(x1, y1, 'r', label='Trajectoire masse 1')
    plt.plot(x2, y2, 'g', label='Trajectoire masse 2')
    plt.xlabel('Position X')
    plt.ylabel('Position Y')
    plt.title(f'Trajectoire du pendule pour alpha = {alpha}')
    plt.legend()
    plt.axis('equal')
    plt.savefig(filename)
    plt.close()

# Fonction pour tracer les angles theta1 et theta2 en fonction du temps
def angle_plot(temps, sol, alpha, filename):
    theta1 = sol[0, :]
    theta2 = sol[1, :]
    plt.figure()
    plt.plot(temps, np.degrees(theta1), label=r'$\theta_1(t)$')
    plt.plot(temps, np.degrees(theta2), label=r'$\theta_2(t)$')
    plt.xlabel('Temps [s]')
    plt.ylabel('Angle [degrés]')
    plt.title(f'Angles des pendules pour alpha = {alpha}')
    plt.legend()
    plt.savefig(filename.format(alpha=alpha))
    plt.close()


# Fonction pour tracer les vitesses angulaires en fonction du temps
def angular_velocity_plot(temps, sol, alpha, filename):
    dth1 = sol[2, :]
    dth2 = sol[3, :]
    plt.figure()
    plt.plot(temps, dth1, label=r'$\dot{\theta}_1(t)$')
    plt.plot(temps, dth2, label=r'$\dot{\theta}_2(t)$')
    plt.xlabel('Temps [s]')
    plt.ylabel('Vitesse angulaire [rad/s]')
    plt.title(f'Vitesses angulaires pour alpha = {alpha}')
    plt.legend()
    plt.savefig(filename.format(alpha=alpha))
    plt.close()


def PhaseSpacePlot(temps, sol, l1, l2, m1, m2, alpha, filename):
    theta1 = sol[0, :]
    dtheta1 = sol[2, :]
    theta2 = sol[1, :]
    dtheta2 = sol[3, :]
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.plot(theta1, dtheta1, 'b', label=r'Espace de phase $\theta_1, \dot{\theta}_1$')
    plt.title(r'Espace de phase pour $\theta_1$ pour alpha = {alpha}')
    plt.xlabel(r'$\theta_1$ (rad)')
    plt.ylabel(r'$\dot{\theta}_1$ (rad/s)')
    plt.grid(True)
    plt.legend()
    plt.subplot(1, 2, 2)
    plt.plot(theta2, dtheta2, 'g', label=r'Espace de phase $\theta_2, \dot{\theta}_2$')
    plt.title(r'Espace de phase pour $\theta_2$ pour alpha = {alpha}')
    plt.xlabel(r'$\theta_2$ (rad)')
    plt.ylabel(r'$\dot{\theta}_2$ (rad/s)')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename.format(alpha=alpha))
    plt.close()


# Fonction pour calculer l'énergie à chaque pas de temps
def energyCompute(sol, l1, l2, m1, m2):
    g = 9.81  # Constante gravitationnelle
    th1, th2, dth1, dth2 = sol
    Ec = (0.5 * m1 * (l1 * dth1) ** 2 + 0.5 * m2 * ((l1 * dth1) ** 2 + (l2 * dth2) ** 2 + 2 * l1 * l2 * dth1 * dth2 * np.cos(th1 - th2)))
    Ep = (-m1 * g * l1 * np.cos(th1) - m2 * g * (l1 * np.cos(th1) + l2 * np.cos(th2)))
    Et = Ec + Ep
    return Ec, Ep, Et

# Fonction pour tracer l'évolution de l'énergie
def energy_evolution_plot(temps, Ec_values, Ep_values, Et_values, alpha, filename):
    plt.figure(figsize=(10, 6))
    plt.plot(temps, Ec_values, label='Energie cinétique $E_c(t)$')
    plt.plot(temps, Ep_values, label='Energie potentielle $E_p(t)$')
    plt.plot(temps, Et_values, label='Energie mécanique $E_m(t)$', linestyle='--')
    plt.xlabel('T (s)')
    plt.ylabel("Energie (Joules)")
    plt.title(f"Evolution de l'énergie mécanique pour alpha = {alpha}")
    plt.legend()
    plt.grid(True)
    plt.savefig(filename)
    plt.close()


# Créer et enregistrer l'animation combinée
create_combined_animation(alpha_values, l1, l2, m1, m2, dt, T, 'combined_pendulums.mp4')

# Simuler pour différentes valeurs d'alpha
for alpha in alpha_values:
    Y0 = [alpha, alpha, 0, 0]
    temps, sol = method_RK4(Y0, l1, l2, m1, m2, dt, T)
    print(f"\nProcessed alpha = {alpha}")

    # Calcul des énergies
    energies = np.array([energyCompute(sol[:, i], l1, l2, m1, m2) for i in range(sol.shape[1])])
    Ec_values = energies[:, 0]
    Ep_values = energies[:, 1]
    Et_values = energies[:, 2]

    # Enregistrement des courbes d'énergie
    energy_evolution_plot(temps, Ec_values, Ep_values, Et_values, alpha, f'energy_evolution_alpha_{alpha}.png')

    # Vérification de la conservation de l'énergie
    energy_conserved = np.std(Et_values) < 1e-5
    energy_conservation_status = 'conservée' if energy_conserved else 'non conservée'
    print(f"L'énergie du double pendule pour alpha = {alpha} est {energy_conservation_status}.")
    print(f"La différence d'énergie mécanique est {np.std(Et_values)} pour alpha = {alpha}")

    plot_trajectories(temps, sol, f'trajectories_alpha_{alpha}.png', alpha)
    angle_plot(temps, sol, alpha, 'angle_plot_alpha_{alpha}.png')
    angular_velocity_plot(temps, sol, alpha, 'angular_velocity_plot_alpha_{alpha}.png')
    PhaseSpacePlot(temps, sol, l1, l2, m1, m2, alpha, 'phase_space_plot_alpha_{alpha}.png')
