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

def test():
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

def qgrid_level_sets():
    """Create a grid of
    (q_1, q_2) plots showing
    level sets of H for
    different values of p
    (the rows/cols of the
    plots vary the p value)
    """
    N = 100
    n_plots = 3
    mu = 0.25
    qmin = -2**0.5
    qmax = -qmin
    pmin = -1
    pmax = 1

    q1s = np.linspace(qmin, qmax, N)
    q2s = np.linspace(qmin, qmax, N)
    p1s = np.linspace(pmin, pmax, n_plots)
    p2s = np.linspace(pmin, pmax, n_plots)

    Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s)
    H = Hamiltonian(Q1, Q2, P1, P2, mu=mu)

    fig, axs = plt.subplots(n_plots, n_plots)
    fig.suptitle("Level sets of the Hamiltonian $H_{%s}$" % mu)
    for i in range(n_plots):
        for j in range(n_plots):
            axs[i,j].contour(
                Q1[:,:,i,j],
                Q2[:,:,i,j],
                H[:,:,i,j],
                levels=np.linspace(-2, 0, 20)
            )
            axs[i,j].set_title("$p = (%s, %s)$" % (P1[0,0,i,j], P2[0,0,i,j]))

def pgrid_level_sets():
    """Create a grid of
    (p_1, p_2) plots showing
    level sets of H for
    different values of p
    (the rows/cols of the
    plots vary the p value)
    """
    N = 100
    n_plots = 5
    mu = 0.25
    qmin = -2**0.5
    qmax = -qmin
    pmin = -1
    pmax = 1

    q1s = np.linspace(qmin, qmax, n_plots)
    q2s = np.linspace(qmin, qmax, n_plots)
    p1s = np.linspace(pmin, pmax, N)
    p2s = np.linspace(pmin, pmax, N)

    Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s)
    H = Hamiltonian(Q1, Q2, P1, P2, mu=mu)

    fig, axs = plt.subplots(n_plots, n_plots)
    fig.suptitle("Level sets of the Hamiltonian $H_{%s}$" % mu)
    for i in range(n_plots):
        for j in range(n_plots):
            axs[i,j].contour(
                P1[i,j,:,:],
                P2[i,j,:,:],
                H[i,j,:,:],
                levels=np.linspace(-2, 0, 20)
            )
            axs[i,j].set_title("$q = (%s, %s)$" % (Q1[i,j,0,0], Q2[i,j,0,0]))


def side_view_q1p1():
    N = 100
    n_plots = 1
    mu = 0.25
    qmin = -2**0.5
    qmax = -qmin
    pmin = -2
    pmax = 2

    q1s = np.linspace(qmin, qmax, N)
    #q2s = np.linspace(qmin, qmax, n_plots)
    q2s = np.array([0])
    p1s = np.linspace(pmin, pmax, N)
    #p2s = np.linspace(pmin, pmax, N)
    p2s = np.array([0])

    print(q1s.shape, q2s.shape, p1s.shape, p2s.shape)
    Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s, indexing="ij")
    H = Hamiltonian(Q1, Q2, P1, P2, mu=mu)

    fig, ax = plt.subplots()
    ax.contour(
                Q1[:,0,:,0],
                P1[:,0,:,0],
                H[:,0,:,0],
                levels=np.linspace(-2, 0, 20)
            )

def side_views_q1p1_p2_slices():
    N = 100
    n_plots = 5
    mu = 0.25
    qmin = -2**0.5
    qmax = -qmin
    q2min = 0
    q2max = 0.2
    pmin = -2
    pmax = 2

    q1s = np.linspace(qmin, qmax, N)
    q2s = np.linspace(q2min, q2max, int(n_plots**0.5)**2)
    p1s = np.linspace(pmin, pmax, N)
    #p2s = np.linspace(pmin, pmax, N)
    p2s = np.array([0])

    Q1, Q2, P1, P2 = np.meshgrid(q1s, q2s, p1s, p2s, indexing="ij")
    H = Hamiltonian(Q1, Q2, P1, P2, mu=mu)

    fig, axs = plt.subplots(*2*(int(n_plots**0.5), ), squeeze=True)
    fig.suptitle("View in the $q_1$-$p_1$ plane of level sets of the Hamiltonian, with $p_2 = 0$ and various values of $q_2$")
    for i, ax in enumerate(np.ravel(axs)):
        ax.set_title("$q_2 = %s$" % Q2[0,i,0,0])
        ax.set_xlabel("$q_1$")
        ax.set_ylabel("$p_1$")
        ax.contour(
                    Q1[:,i,:,0],
                    P1[:,i,:,0],
                    H[:,i,:,0],
                    levels=np.linspace(-2, 0, 20)
                )

