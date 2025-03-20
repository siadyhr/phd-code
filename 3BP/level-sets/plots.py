import numpy as np
import matplotlib.pyplot as plt

def V(q1, q2, mu=0.5):
    return -(
            mu/np.sqrt(
                (q1 + mu)**2 + q2**2
            )
            +
            (1-mu)/np.sqrt(
                (q1 + mu - 1)**2 + q2**2
            )
        )

def Hamiltonian(q1, q2, p1, p2, mu=0.5):
    return (
            (p1**2 + p2**2)/2
            +
            q1*p2 - q2*p1
            +
            V(q1, q2, mu=mu)
        )

N = 100
q1s = np.linspace(-np.sqrt(2), np.sqrt(2), N)
q2s = np.linspace(-np.sqrt(2), np.sqrt(2), N)
p1s = np.array([-0.2])
p2s = np.array([0])

Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s)
print(Q1.shape)

H = Hamiltonian(Q1, Q2, P1, P2, mu=0.25)

fig, ax = plt.subplots()
#ax.plot(qs[:,0], V(qs)) 
print(np.min(H), np.max(H))
ax.contour(Q1[:,:,0,0], Q2[:,:,0,0], H[:,:,0,0], levels=np.linspace(-5, 0, 50))

# Idé 1: p-gitter af subplots med level sets
# Idé 2: Level-set for fast energi, men tegn forskellige med forskellig p (eller mu!)

# Idé 1
N = 100
n_plots = 3
q1s = np.linspace(-np.sqrt(2), np.sqrt(2), N)
q2s = np.linspace(-np.sqrt(2), np.sqrt(2), N)
p1s = np.linspace(-1, 1, n_plots)
p2s = np.linspace(-1, 1, n_plots)

Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s)
H = Hamiltonian(Q1, Q2, P1, P2, mu=0.25)

fig, axs = plt.subplots(n_plots, n_plots)
for i in range(n_plots):
    for j in range(n_plots):
        axs[i,j].contour(
            Q1[:,:,i,j],
            Q2[:,:,i,j],
            H[:,:,i,j],
            levels=np.linspace(-2, 0, 20)
        )

plt.show()
