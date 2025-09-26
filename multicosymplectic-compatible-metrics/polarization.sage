def is_compatible(alpha, beta, g):
    assert (beta + beta.T == 0), "beta is not antisymmetric"
    assert (beta.kernel().dimension() == 1), "dim(ker beta) ≠ 1"
    R = beta.kernel().basis()[0]    # is vector
    R = matrix(R).T # to get R as column
    R /= (alpha * R)
    phi_candidate = g^(-1) * beta
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

