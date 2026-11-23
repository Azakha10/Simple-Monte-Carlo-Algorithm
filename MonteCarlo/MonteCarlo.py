import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Q1. Linear Congruential Generator #
def lcg(seed,a,c,m,n):
  liste = []
  for i in range(n):
    seed = (a*seed+c) % m
    liste.append(seed/m)
  return liste

# Q2. MAIN #
liste = lcg(42,1664525,1013904223,2**32,10000)
print(liste)
print("Longueur :", len(liste))

# Q3. X2 #
classe1 = []
classe2 = []
classe3 = []
classe4 = []
classe5 = []

for x in liste:
  if 0 <= x < 0.2:
    classe1.append(x)
  elif 0.2 <= x < 0.4:
    classe2.append(x)
  elif 0.4 <= x < 0.6:
    classe3.append(x)
  elif 0.6 <= x < 0.8:
    classe4.append(x)
  elif 0.8 <= x < 1.0:
    classe5.append(x)

# Q1. Proba d'occurence #
sigma_photo = 4.37
sigma_rayleigh = 1.30
sigma_compton = 3.19
sigma_paires = 0.0
sigma_total = sigma_photo + sigma_rayleigh + sigma_compton + sigma_paires

p_photo = sigma_photo / sigma_total
p_rayleigh = sigma_rayleigh / sigma_total
p_compton = sigma_compton / sigma_total
p_paires = sigma_paires / sigma_total

# Q1. Methode Directe (Distance) #
def tirage_distance(lambda_mean):
  s = []
  for i in range(len(liste)):
    s.append(-lambda_mean*np.log(1-liste[i]))
  return s

# Q2. Methode Directe (Direction) #
def tirage_direction():
  theta = []
  phi = []
  for i in range(len(liste)):
    theta.append(np.arccos(1 - 2*liste[i]))
    phi.append(2*np.pi*liste[i])
  return theta,phi

# Q3. Methode de Rejet #
def tirage_rejet():
  echantillons_acceptes = []
  essais = len(liste) // 2

  for i in range(essais):
    R1 = liste[2 * i]
    xi = np.pi * R1
    R2 = liste[2 * i + 1]
    if xi < 1e-9: # Pour éviter division par zero
        acceptance = 1.0
    else:
        acceptance = np.sin(xi) / xi
    if R2 <= acceptance:
      echantillons_acceptes.append(xi)
  return echantillons_acceptes

# Q2. Choix d'interaction #
def type_interaction(r):
    if r < p_photo:
      return "photo"
    elif r < (p_photo + p_rayleigh):
      return "reyl"
    elif r < (p_photo + p_rayleigh + p_compton):
      return "comp"
    else:
      return "pair"

# Q3. Algorithme principal 3D #
def transport_photon_3d(lambda_mean):
    distance_totale = 0
    nombre_interactions = 0
    absorbe = False
    i = 0

    # Listes pour stocker la trajectoire 3D
    pos_x, pos_y, pos_z = [0.0], [0.0], [0.0]
    liste2 = lcg(seed, 1664525, 1013904223, 2**32, 40000)

    while not absorbe and nombre_interactions < 10000:
        r_distance = liste2[i]
        r_theta = liste2[i + 1]
        r_phi = liste2[i + 2]
        r_interaction = liste2[i + 3]
        i += 4
        nombre_interactions += 1

        # 1. Tirage de la distance et des angles
        d = -lambda_mean * np.log(1 - r_distance)
        distance_totale += d

        # Tirage theta et phi pour la direction 3D
        theta = np.arccos(1 - 2 * r_theta)
        phi = 2 * np.pi * r_phi

        # Conversion coordonnées sphériques en cartésiennes
        dx = d * np.sin(theta) * np.cos(phi)
        dy = d * np.sin(theta) * np.sin(phi)
        dz = d * np.cos(theta)

        # Mise à jour des coordonnées
        pos_x.append(pos_x[-1] + dx)
        pos_y.append(pos_y[-1] + dy)
        pos_z.append(pos_z[-1] + dz)

        # 2. Choix du type d'interaction
        interaction = type_interaction(r_interaction)

        if interaction == "photo":
            absorbe = True
            print(f"Le photon est absorbé après {nombre_interactions} interactions. (Effet Photoelectrique)")
            print(f"Distance totale parcourue : {distance_totale:.4f} m")
        else:
            print(f"{nombre_interactions}eme interaction : {interaction}")

    # Tracé de la trajectoire en 3D
    print("")
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(pos_x, pos_y, pos_z, 'k-o', label='Trajectoire du photon')
    ax.scatter(pos_x[0], pos_y[0], pos_z[0], color='green', s=100, label='Source')
    ax.scatter(pos_x[-1], pos_y[-1], pos_z[-1], color='red', marker='x', s=100, label='Absorption')

    ax.set_title('Trajectoire 3D du photon')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_zlabel('Z (m)')
    ax.legend()
    plt.show()

    return distance_totale

seed = 1000
print("---------------------------------")
print("--- Simulation 3D d'un photon ---")
print("---------------------------------")
transport_photon_3d(1)
print(f"(Seed : {seed})")