import matplotlib
import tkinter as tk

matplotlib.use('TkAgg')  # Utilisez un backend interactif
import os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

plt.rc('figure', figsize=(12, 9))
plt.rcParams['animation.writer'] = 'ffmpeg'


def F(Y, l1, l2, m1, m2):
    th1, th2, dth1, dth2 = Y
    g = 9.81  # Constante gravitationnelle

    d2th1 = (-dth1**2*l1*m2*np.sin(2*th1 - 2*th2)/2 - dth2**2*l2*m2*np.sin(th1 - th2) - g*m1*np.sin(th1) - g*m2*np.sin(th1)/2 - g*m2*np.sin(th1 - 2*th2)/2)/(l1*(m1 - m2*np.cos(th1 - th2)**2 + m2))
    d2th2 = (dth1**2*l1*m1*np.sin(th1 - th2) + dth1**2*l1*m2*np.sin(th1 - th2) + dth2**2*l2*m2*np.sin(2*th1 - 2*th2)/2 - g*m1*np.sin(th2)/2 + g*m1*np.sin(2*th1 - th2)/2 - g*m2*np.sin(th2)/2 + g*m2*np.sin(2*th1 - th2)/2)/(l2*(m1 - m2*np.cos(th1 - th2)**2 + m2))

    return np.array([dth1, dth2, d2th1, d2th2])


#Euler explicite
def method_Euler_explicite(Y, l1, l2, m1, m2, dt, T):
    N = int(T / dt)
    temps = np.linspace(0, T, N)
    sol = np.zeros((4, N))
    sol[:, 0] = Y
    for i in range(N - 1):
        sol[:, i + 1] = sol[:, i] + dt * F(sol[:, i], l1, l2, m1, m2)
    return temps, sol


# Affichage
def visualisation(temps, sol, speed_up_factor=1000):
    fig, ax = plt.subplots()
    x1 = l1 * np.sin(sol[0, :])
    y1 = -l1 * np.cos(sol[0, :])
    x2 = x1 + l2 * np.sin(sol[1, :])
    y2 = y1 - l2 * np.cos(sol[1, :])
    ax.set_xlim(-(l1 + l2), l1 + l2)
    ax.set_ylim(-(l1 + l2), l1 + l2)

    # Initialisation de la ligne pour les barres du pendule
    line, = plt.plot([], [], 'o-', linewidth=2, color="blue")
    # Initialisation des trajectoires pour chaque masse
    path1, = plt.plot([], [], ',', color="red")  # Trajectoire pour la première masse
    path2, = plt.plot([], [], ',', color="green")  # Trajectoire pour la deuxième masse

    # Stockage des positions pour les trajectoires
    x1data, y1data = [], []
    x2data, y2data = [], []

    def maj(i):
        # Mise à jour des barres du pendule
        line.set_data([0, x1[i], x2[i]], [0, y1[i], y2[i]])
        # Ajout des positions actuelles aux trajectoires
        x1data.append(x1[i])
        y1data.append(y1[i])
        x2data.append(x2[i])
        y2data.append(y2[i])
        # Mise à jour des trajectoires
        path1.set_data(x1data, y1data)
        path2.set_data(x2data, y2data)
        print(f"Enregistrement de l'image {i + 1}/{len(temps)}")
        return line, path1, path2,

    interval_ms = 50 / speed_up_factor
    ani = animation.FuncAnimation(fig, maj, frames=range(len(x1)), blit=True, interval=interval_ms)
    return ani

Y0 = [np.pi/6, np.pi/6, 0, 0]  # [th1, th2, dth1, dth2]

# Constantes du système
m1, m2, l1, l2 = 1, 1, 1, 1

# Application
N = 1000
temps_aquis = 1
dt = temps_aquis / N

temps, sol = method_Euler_explicite(Y0, l1, l2, m1, m2, dt, T=temps_aquis)
ani = visualisation(temps, sol)


# Fonction pour tracer les angles theta1 et theta2 en fonction du temps
def angle_plot(temps, sol, filename='angle_plot.png'):
    theta1 = sol[0, :]
    theta2 = sol[1, :]

    plt.figure()
    plt.plot(temps, np.degrees(theta1), label=r'$\theta_1(t)$')
    plt.plot(temps, np.degrees(theta2), label=r'$\theta_2(t)$')
    plt.xlabel('Temps [s]')
    plt.ylabel('Angle [degrés]')
    plt.title('Angles des pendules en fonction du temps')
    plt.legend()
    plt.savefig(filename)
    plt.close()


# Fonction pour tracer les vitesses angulaires en fonction du temps
def angular_velocity_plot(temps, sol, filename='angular_velocity_plot.png'):
    dth1 = sol[2, :]
    dth2 = sol[3, :]

    plt.figure()
    plt.plot(temps, dth1, label=r'$\dot{\theta}_1(t)$')
    plt.plot(temps, dth2, label=r'$\dot{\theta}_2(t)$')
    plt.xlabel('Temps [s]')
    plt.ylabel('Vitesse angulaire [rad/s]')
    plt.title('Vitesses angulaires des pendules en fonction du temps')
    plt.legend()
    plt.savefig(filename)
    plt.close()


# Fonction pour tracer et enregistrer les trajectoires des masses du pendule
def plot_trajectories(temps, sol, filename='trajectories.png'):
    # Calcul des positions des masses
    theta1 = sol[0, :]
    theta2 = sol[1, :]
    x1 = l1 * np.sin(theta1)
    y1 = -l1 * np.cos(theta1)
    x2 = x1 + l2 * np.sin(theta2)
    y2 = y1 - l2 * np.cos(theta2)
    plt.figure()
    plt.plot(x1, y1, 'r', label='Trajectoire de la première masse')
    plt.plot(x2, y2, 'g', label='Trajectoire de la deuxième masse')
    plt.xlabel('Position X')
    plt.ylabel('Position Y')
    plt.title('Trajectoires des masses du pendule')
    plt.legend()
    plt.axis('equal')  # Assurer que les proportions sont égales
    plt.savefig(filename)
    plt.close()


# Chemin absolu du script courant
chemin_script_courant = os.path.abspath(__file__)

# Répertoire contenant le script courant
repertoire_script_courant = os.path.dirname(os.path.abspath(__file__))

# Application de la visualisation et des tracés
#angle_plot(temps, sol, filename=os.path.join(repertoire_script_courant, 'angle_plot_2.png'))
#angular_velocity_plot(temps, sol, filename=os.path.join(repertoire_script_courant, 'angular_velocity_plot_2.png'))
#plot_trajectories(temps, sol, filename=os.path.join(repertoire_script_courant, 'trajectories_2.png'))

def save_animation_once(ani, filename='double_pendulum_2sec_2.mp4', fps=200):
    if not hasattr(save_animation_once, "already_saved"):  # Vérifie si la fonction a déjà été appelée
        ani.save(filename, fps=fps)
        print(f"Animation enregistrée sous : {filename}")
        save_animation_once.already_saved = True  # Marque la fonction comme déjà appelée

# Pour afficher l'animation dans une fenêtre (sans enregistrer)
plt.show()

# Pour sauvegarder l'animation
#save_animation_once(ani, os.path.join(repertoire_script_courant, 'double_pendulum_2sec_2.mp4'))

def energyCompute(sol, l1, l2, m1, m2):
    g = 9.81  # Constante gravitationnelle
    th1, th2, dth1, dth2 = sol
    #Energie cinétique
    Ec = (0.5 * m1 * (l1 * dth1) ** 2 + 0.5 * m2 * ((l1 * dth1) ** 2 + (l2 * dth2) ** 2 + 2 * l1 * l2 * dth1 * dth2 * np.cos(th1 - th2)))
    #Eenrgie potentielle
    Ep = (-m1 * g * l1 * np.cos(th1) - m2 * g * (l1 * np.cos(th1) + l2 * np.cos(th2)))
    #Energie mécanique
    Et = Ec + Ep
    return Ec, Ep, Et

# Calcule les énergies
energies = np.array([energyCompute(sol[:, i], l1, l2, m1, m2) for i in range(sol.shape[1])])
Ec_values = energies[:, 0]
Ep_values = energies[:, 1]
Et_values = energies[:, 2]

def energy_evolution_plot(temps, Ec_values, Ep_values, Et_values, filename):
    plt.figure(figsize=(10, 6))
    plt.plot(temps, Ec_values, label='Energie cinétique $E_c(t)$')
    plt.plot(temps, Ep_values, label='Eenrgie potentielle $E_p(t)$')
    plt.plot(temps, Et_values, label='Energie mécanique $E_m(t)$', linestyle='--')
    plt.xlabel('T(s)')
    plt.ylabel('Energie (Joules)')
    plt.title("Evolution de l'énergie mécanique du double pendule")
    plt.legend()
    plt.grid(True)
    plt.savefig(filename)
    plt.close()

energy_evolution_plot(temps, Ec_values, Ep_values, Et_values, filename=os.path.join(repertoire_script_courant, 'energy_evolution_2.png'))

energy_conserved = np.std(Et_values) < 1e-5

# Indicate whether the energy is conserved
energy_conservation_status = 'conservée' if energy_conserved else 'non conservée'
print(f"L'énergie du double pendule est {energy_conservation_status}.")
print(f"La différence d'énergie mécanique est {np.std(Et_values)}")
