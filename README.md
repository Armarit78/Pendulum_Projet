# Modélisation du Chaos du Pendule Double

> **Note :** Pour plus d'informations détaillées sur le projet, veuillez consulter le rapport complet disponible dans ce dépôt.

## Aperçu du Projet
Ce projet explore le monde fascinant des systèmes chaotiques à travers l'exemple du pendule double. Connus pour leur comportement imprévisible et complexe, les pendules doubles illustrent parfaitement comment de légères variations des conditions initiales peuvent entraîner des résultats radicalement différents. Cette sensibilité fait du pendule double un sujet d'étude classique dans les domaines de la théorie du chaos et des systèmes non linéaires.

## Modèle Mathématique
Le système est constitué de deux pendules attachés bout à bout, où le mouvement est régi par les lois de la mécanique classique. La position des masses est décrite par les angles \( \theta_1 \) et \( \theta_2 \), qui déterminent les trajectoires des pendules sous l'effet de la gravité.

## Méthodes Numériques
Nous avons abordé le problème en utilisant des méthodes numériques pour résoudre les équations différentielles ordinaires dérivées du modèle physique :
- **Méthode d'Euler :** Une méthode directe pour l'intégration numérique.
- **Méthode de Runge-Kutta (Ordre 4) :** Offre une précision et une stabilité supérieures.
- **Intégration de Verlet :** Idéale pour les systèmes où la conservation de l'énergie est critique.

## Contenu du Dépôt
- **Partie 1 :** Contient les implémentations et exemples utilisant la méthode explicite d'Euler.
- **Partie 2 :** Contient les implémentations et exemples utilisant la méthode de Verlet.
- **Partie 3 :** Contient des analyses sur le mouvement chaotique et la sensibilité aux conditions initiales.

## Contributeurs
- Guillaume PORET
- Ludovic CURE-MOOG
