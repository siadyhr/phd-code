# Standard form
# alpha = dx^5
# kerne er <d/dx^1, ..., d/dx^4>
beta1 = matrix([
    [ 0, 1, 0, 0, 0], 
    [-1, 0, 0, 0, 0], 
    [ 0, 0, 0, 1, 0], 
    [ 0, 0,-1, 0, 0], 
    [ 0, 0, 0, 0, 0]
])
alpha1 = matrix([0,0,0,0,1])

beta2 = matrix([
    [ 0, 0, 0, 0, 0],
    [ 0, 0, 1, 0, 0],
    [ 0,-1, 0, 0, 0],
    [ 0, 0, 0, 0, 1],
    [ 0, 0, 0,-1, 0],
])
alpha2 = matrix([1,1,0,0,0])

G0 = diagonal_matrix([2,1,1,1,1])

def Reeb_field(alpha, beta):
    R = beta.kernel().basis()[0]    # is vector
    R = matrix(R).T # to get R as column
    R /= (alpha * R)
    return R

def is_compatible(alpha, beta, g):
    assert (beta + beta.T == 0), "beta is not antisymmetric"
    assert (beta.kernel().dimension() == 1), "dim(ker beta) ≠ 1"
    phi_candidate = g^(-1) * beta
    R = Reeb_field(alpha, beta)
#    print("φ² = ")
#    print(phi_candidate**2)
#    print("|R|^2 =",(R.T * g * R)[0,0])
    return (
            (
                phi_candidate**2
                ==
                -identity_matrix(beta.dimensions()[0])
                +
                R * alpha
            ) and (
                R.T * g * R == 1
            )
    )

def diagonalize(A):
    eigendata = A.eigenvectors_right()
    eigenvectors = []
    eigenvalues = []
    for (eigenvalue, eigenvector_list, multiplicity) in eigendata:
        for _ in range(multiplicity):
            eigenvalues.append(eigenvalue)
        for eigenvector in eigenvector_list:
            eigenvectors.append(eigenvector)
    D = diagonal_matrix(eigenvalues)
    V = matrix(eigenvectors)
    return (D, V)

def inverse_square_root(A):
    D, V = diagonalize(A)
    Dsqrt = diagonal_matrix(
            x**(-1/2) if x != 0 else 0
            for x in D.diagonal()
    )
    return (1/V) * Dsqrt * V

def polarize(alpha, beta, g0):
    """
    Mathematical procedure:
        beta = g0(-, A -)
    """

    A = g0^(-1) * beta
    J = (inverse_square_root(A * A.T)) * A
    g = J.T * beta + alpha.T * alpha
    return g

def print_structure(alpha, beta):
    print("Cosymplectic structure:")
    print("α = %s" % alpha)
    print("β =")
    print(beta)
    print("R = %s" % Reeb_field(alpha, beta).T)

print("Structure 1")
print_structure(alpha1, beta1)
print("Structure 2")
print_structure(alpha2, beta2)

print("is α_i(R_j) = δ_ij?")
print(alpha2 * Reeb_field(alpha1, beta1))
print(alpha1 * Reeb_field(alpha2, beta2))

print("Is G₀ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G0))
G1 = polarize(alpha1, beta1, G0)
print(G1)
print("Is G₁ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G1))

print()
print("Polarize G₁ for (α₂, β₂). g₂ =")
G2 = polarize(alpha2, beta2, G1)
print(G2)
print("Is G₂ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G2))
print("Is G₂ compatible with (α₂, β₂)?", is_compatible(alpha2, beta2, G2))
