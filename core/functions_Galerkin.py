import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Galerkin:
    """
    Class TP2 Analyse numérique
    """
    def __init__(self, L, n):
        """
        Initialisation de la classe Galerkin.

        Paramètres :
        - L : Longueur totale du domaine
        - n : Nombre de points de discrétisation
        """
        self.L = L
        self.n = n
        self.h = L / n
        self.sous_diag, self.diag, self.sur_diag, self.S, self.s_2 = self.initialize_matrices()
        self.c_tilde_trapeze=self.algo_thomas(self.S)
        self.c_tilde_exate=self.algo_thomas(self.s_2)


    def initialize_matrices(self):
        """
        Initialisation des matrices pour le système linéaire tridiagonal.

        Sortie :
        - sous_diag, diag, sur_diag : Les sous-diagonales, diagonales et sur-diagonales de la matrice tridiagonale du système linéaire
        - S : second membre du système linéaire
        """
        sous_diag = np.zeros(self.n + 1)
        diag = np.zeros(self.n + 1)
        sur_diag = np.zeros(self.n + 1)
        S = np.zeros(self.n + 1)

        i_values = np.arange(1, self.n + 1)
        sous_diag[1:self.n+1] = (-1 / self.h) - (i_values - 1) + (1 / 2)
        diag[1:self.n+1] = (2 / self.h) + 2 * (i_values - 1)
        sur_diag[1:self.n+1] = (-1 / self.h) - (i_values - 1) - (1 / 2)

        sous_diag[0] = 0
        sous_diag[self.n] = 1/self.h + - self.n + (1 / 2)

        diag[0] = (1 / self.h) + (1 / 2)
        diag[self.n] = 1 / self.h + self.n - 1 / 2

        sur_diag[0] = 1 / (2 * self.h)
        sur_diag[self.n] = 1/2


        # Second membre avec la méthode des trapèzes

        x = np.arange(0, self.L, self.h)
        S[1:self.n+1]=self.h*(x)**2
        S[0] = 0
        S[self.n] = self.h/2

        # Second membre avec le calcul exate de l'intégrale
        s_2 = np.zeros(self.n + 1)
        s_2[1:self.n+1] = (2 * ((i_values - 1) * self.h)**4 - (self.h * i_values)**4 - (self.h * (i_values - 2))**4) * (1 / (4 * self.h) - 1 / (3 * self.h))
        s_2[0] = (self.h**3) / 12
        s_2[self.n] = (1 / self.h) * ((self.L**4) / 4 - (self.L**3) / 3 * (self.n - 1) * self.h + (((self.n - 1) * self.h)**4) / 12)

        return sous_diag, diag, sur_diag, S, s_2

    def algo_thomas(self, s):
        """
        Algorithme de Thomas pour la résolution d'un système linéaire tridiagonal.
        Entrée : 
        - s : second membre du système linéaire
        Sortie :
        - c_tilde : Solution du système linéaire
        """
        c_thomas = np.zeros(self.n + 1)
        d_thomas = np.zeros(self.n + 1)
        c_tilde = np.zeros(self.n + 1)

        c_thomas[0] = self.sur_diag[0] / self.diag[0]
        for i in range(1, self.n):
            c_thomas[i] = self.sur_diag[i] / (self.diag[i] - c_thomas[i - 1] * self.sous_diag[i])

        d_thomas[0] = s[0] / self.diag[0]
        for i in range(1, self.n + 1):
            d_thomas[i] = (s[i] - d_thomas[i - 1] * self.sous_diag[i]) / (self.diag[i] - c_thomas[i - 1] * self.sous_diag[i])

        c_tilde[self.n] = d_thomas[self.n]
        for i in range(1, self.n):
            c_tilde[self.n - i] = d_thomas[self.n - 1 - i] - c_thomas[self.n - 1 - i] * c_tilde[self.n - i + 1]

        c_tilde[0] = 0
        c_tilde[self.n] = 0
        return c_tilde

    def solution_exate(self, x):
        """
        Calcul de la solution exacte.

        Paramètres :
        - x : Valeurs de x où la solution exacte doit être évaluée

        Sortie :
        - Résultat de la solution exacte pour chaque valeur de x
        """
        return -1/3 * x + 1/6 * x**2 - 1/9 * x**3 + 0.4007486225 * np.log(1 + x)


    
    
    def plot_results(self):
        """
        Tracé des résultats obtenus par la méthode de Galerkin et la solution exacte.
        """
        sns.set(style="whitegrid")
        x = np.arange(0, self.L + self.h, self.h)
        c_exact = self.solution_exate(x)
        
        plt.plot(x, self.c_tilde_trapeze, 'r', label='$c_{num}$')
        plt.plot(x, c_exact, label="$c_{ex}$")
        plt.legend()
        plt.show()
    
    def calculate_error(self):
        """
        Calcul de l'erreur quadratique entre la solution obtenue par la méthode de Galerkin
        et la solution exacte.

        Sortie :
        - Erreur quadratique
        """

        x = np.arange(0, self.L + self.h, self.h)
        c_exact = self.solution_exate(x)
        error_trapeze=np.sqrt(np.mean((c_exact - self.c_tilde_trapeze)**2))
        error_exate=np.sqrt(np.mean((c_exact - self.c_tilde_exate)**2))
        return error_trapeze, error_exate

 

