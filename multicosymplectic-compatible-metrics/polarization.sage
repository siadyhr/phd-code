import numpy as np
import scipy
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

#G0 = diagonal_matrix([2,1,1,1,1])
G0 = random_matrix(QQ, 5)

def Reeb_field(alpha, beta):
    R = beta.kernel().basis()[0]    # is vector
    R = matrix(R).T # to get R as column
    R /= (alpha * R)
    return R

def is_cosymplectic(alpha, beta):
    return (
            (beta == -beta.T)
            and
            (beta.nullity() == 1)
            and
            (
                (alpha * (beta.kernel().basis()[0]))[0]
                !=
                0
            )
            )

def is_compatible(alpha, beta, g):
    assert (beta + beta.T == 0), "beta is not antisymmetric"
    assert (beta.kernel().dimension() == 1), "dim(ker beta) ≠ 1"
    phi_candidate = g^(-1) * beta
    R = Reeb_field(alpha, beta)
    print("Checking compatibility...")
    print("φ² - (-1 + R×α) = ")
    print(phi_candidate**2 + identity_matrix(beta.dimensions()[0]) - R * alpha)
    print("|R|^2 =",(R.T * g * R)[0,0])
    print(R * alpha)
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

def diagonalize(A, mode="exact"):
#    print("Begin diagonalization, mode %s" % mode)
    if mode == "numpy":
        eigendata = np.linalg.eig(A)
#        D = diagonal_matrix(np.round(eigendata[0], 16))
        D = diagonal_matrix(eigendata[0])
        V = matrix(eigendata[1].tolist())
#        print("D")
#        print(D)
#        print("V")
#        print(V)
#        print("1/V")
#        print(1/V)
#        print("Residuals")
#        print((1/V) * D * V - A)
    elif mode == "exact":
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

def inverse_square_root(A, mode="exact"):
#    print("Inverse square root? (mode %s)" % mode)
    if mode in ["exact", "numpy"]:
#        print(A.eigenvalues())
        D, V = diagonalize(A, mode)
#        print("Got eigenvalues")
#        print(D.diagonal())
        Dsqrt = diagonal_matrix(
                x**(-1/2) if abs(x) > min(map(abs, D.diagonal())) else 0
                for x in D.diagonal()
        )
#        print("Residuals of square root")
        result = V * Dsqrt * (1/V)
#        print()
        res = result**2*A - identity_matrix(5)
#        print(res)
#        print("Σ =", norm(res))
        return result
    elif mode == "scipy":
#        print("Matrix = ")
#        print(np.array(A))
#        print("Eigenvalues:")
#        print(A.eigenvalues())
        print([x**(0.5) for x in A.eigenvalues()])
        # This cannot work since sqrt(A) should not be invertible...
        # (-> numerical errors)
#        result = np.linalg.inv(scipy.linalg.sqrtm(np.array(A)))
        result = scipy.linalg.sqrtm(scipy.linalg.pinv(A))
        print("Inverse square root = ")
        print(result)
#        print("Eigenvalues of it")
#        print(matrix(result).eigenvalues())
#        print("Residuals =")
        result = matrix(result)
        res = (result**2) * A - identity_matrix(5)
#        print(res)
#        print("Σ =", norm(res))
        return result

def polarize(alpha, beta, g0, mode="exact"):
    """
    Mathematical procedure:
        beta = g0(-, A -)
    """

#    print("Begin polarizing")
#    print("Calculate A")
    A = g0^(-1) * beta
#    print("Calculate J")
    J = (inverse_square_root(A * A.T, mode)) * A
#    print("Calculate g")
    g = J.T * beta + alpha.T * alpha
#    print("Done polarizing")
    return g

def print_structure(alpha, beta):
    print("Cosymplectic structure:")
    print("α = %s" % alpha)
    print("β =")
    print(beta)
    print("R = %s" % Reeb_field(alpha, beta).T)

print("Structure 1")
print_structure(alpha1, beta1)
print("Is cosymplectic?", is_cosymplectic(alpha1, beta1))

print("Structure 2")
print_structure(alpha2, beta2)
print("Is cosymplectic?", is_cosymplectic(alpha2, beta2))

print("is α_i(R_j) = δ_ij?")
print(alpha2 * Reeb_field(alpha1, beta1))
print(alpha1 * Reeb_field(alpha2, beta2))

print("Is G₀ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G0))
G1 = polarize(alpha1, beta1, G0, "numpy")
print(G1)
print("Is G₁ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G1))
quit()

print()
print("Polarize G₁ for (α₂, β₂). g₂ =")
G2 = polarize(alpha2, beta2, G1)
print(G2)
print("Is G₂ compatible with (α₁, β₁)?", is_compatible(alpha1, beta1, G2))
print("Is G₂ compatible with (α₂, β₂)?", is_compatible(alpha2, beta2, G2))
