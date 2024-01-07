# Programme TP1
Programme qui éxécute le TP1

# Installation des bibliothèques Python avec pip

## Prérequis
. Python 3.x installé sur votre machine - https://www.python.org/downloads/


## Installation
1. Ouvrez une invite de commande ou un terminal sur votre machine.


2. Assurez-vous que `pip` est installé en exécutant la commande suivante :

```bash
pip --version
```

Si `pip` n'est pas installé, suivez les étapes de ce tutoriel : 

https://www.youtube.com/watch?v=PikcUT-ts7E&ab_channel=AhmedHegazy 

3. Placez-vous dans le répertoire du projet.

``` bash
cd Chemin_acces
```
Exemple si vous avez téléchargé le fichié sur votre bureau: 

``` bash

cd Desktop\Analyse_numerique_TP1_ailette_de_refroidissement-main\Analyse_numerique_TP1_ailette_de_refroidissement-main
```

4. Installez les bibliothèques nécessaires en exécutant la commande suivante :

``` bash
pip install -r requirements.txt
```
Cette commande va installer toutes les bibliothèques listées dans le fichier requirements.txt qui se trouve dans le répertoire du projet. Assurez-vous que ce fichier contient bien toutes les bibliothèques nécessaires au bon fonctionnement du programme.

5. Lancé le programme


``` bash
python main.py
```