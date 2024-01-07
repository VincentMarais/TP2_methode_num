from core.functions_Galerkin import Galerkin
import matplotlib.pyplot as plt

# Définition des constantes
L = 1
n = int(input("Valeur de n : "))

# Utilisation de la classe Galerkin


galerkin = Galerkin(L, n)

# Affichage des résultats
galerkin.plot_results()

print("Questions prof")
list_n=[2,4,8,65]
error=[]
for i in list_n:
    galerkin = Galerkin(L, i)
    galerkin.plot_results()
    error.append(galerkin.calculate_error()[0])

print("Erreur quadratique question 4:", error)
# Calcul de l'erreur quadratique
error_trapeze = []
error_exact = []
L_values = [2, 4, 20, 65, 100]

for i in L_values:
    galerkin = Galerkin(L, i)
    err_1, err_2 = galerkin.calculate_error()
    error_trapeze.append(err_1)
    error_exact.append(err_2)

plt.plot(L_values, error_trapeze, 'r', label='$c_{trapeze}$')
plt.scatter(L_values, error_trapeze, color='red')
plt.plot(L_values, error_exact, label="$c_{ex}$")
plt.scatter(L_values, error_exact)
plt.legend()
plt.show()


