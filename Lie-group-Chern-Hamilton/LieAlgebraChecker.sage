"""
Functions to check if a collections of tensors
* Form an almost cosymplectic structure
* Is R-invariant if exponentiated from a Lie algebra
* Define a critical metric if exponentiated from a Lie algebra

Regarding 2-tensors and matrices:
    Write S(x, y) as x.T * S * y
    Corresponds to S_{ij} (as a matrix) = S(e_i, e_j)
"""    
debug = 1
def check_compatible_structure(alpha, beta, phi, g):
    """Check if a 1-form, 2-form (1, 1)-tensor and 2-tensor
    actually define an almost cosymplectic structure with a
    compatible metric"""
    assert beta.dimensions() == g.dimensions() == phi.dimensions(), "Inequal dimensions"
    assert len(beta.dimensions()) == 2, "Matrices have dimension %s instead of 2" % beta.dimensions()
    assert beta.dimensions()[0] == beta.dimensions()[1], "Matrices are not square (have dimension %s)" % beta.dimensions()
    assert beta.T + beta == 0, "beta is not antisymmetric"
    assert beta.nullity() == 1, "ker(beta) has dimension %s instead of 1" % beta.nullity()
    R = beta.right_kernel().basis()[0]
    R = matrix(R).T
    assert (alpha * R)[0,0] != 0, "R is in ker(alpha)"
    R = R/((alpha * R)[0,0])
    assert phi**2 == -identity_matrix(beta.dimensions()[0]) + R*alpha, "phi^2 ≠ -1 + alpha ⊗ R"
    assert g*phi == beta, "g(-,phi-) ≠ beta"

def check_alpha_R_invariant(L, alpha, beta, phi, g, R_Liealg):
    """
    R_Liealg
        Element of L representing R
    """
    print("Checking if alpha is R-invariant...")
    for X in L.gens():
        RX = L.bracket(R_Liealg, X).to_vector()
        assert (alpha * RX) == 0
    print("alpha is R-invariant!")

def check_beta_R_invariant(L, alpha, beta, phi, g, R):
    print("Checking if beta is R-invariant...")
    for X in L.gens():
        for Y in L.gens():
            LHS = L.bracket(R, X).to_vector() * beta * Y.to_vector()
            RHS = X.to_vector() * beta * L.bracket(R, Y).to_vector()
            if debug:
                print("Checking (%s, %s)" % (X, Y))
                print("With brackets (%s, %s)" % (L.bracket(R, X), L.bracket(R, Y)))
                print(LHS)
                print(RHS)
            assert (
                    LHS
                    +
                    RHS
                    ==
                    0
                ), "beta is not R-invariant"
    print("beta is R-invariant!")

def check_critical(L, alpha, beta, phi, g, R):
    print("Checking if g is critical...")
    def h(X):
        return (1/2) * (
                L.bracket(R, L.from_vector(phi * X.to_vector()))
                -
                L.from_vector(phi * L.bracket(R, X).to_vector())
            )
    for X in L.gens():
        LHS = L.bracket(
                R, 
                h(L.from_vector(phi * X.to_vector()))
            )
        RHS = h(L.from_vector(phi * L.bracket(R, X).to_vector()))
        assert LHS == RHS, "g is not critical"
    print("g is critical")
    return 1
    # Old (wrong) check
    for X in L.gens():
        for Y in L.gens():
            if debug:
                print("Checking (%s, %s)" % (X, Y))
                print("With brackets (%s, %s)" % (L.bracket(R, X), L.bracket(R, Y)))
                print("phi gives (%s, %s)" % (
                    L.from_vector(phi * X.to_vector()),
                    L.from_vector(phi * Y.to_vector()),
                    ))
            LHS = (
                L.bracket(R, X).to_vector() * g * L.bracket(R, Y).to_vector()
                -
                (
                    L.bracket(
                        R, L.from_vector(phi * X.to_vector())
                    ).to_vector()
                    *
                    g
                    *
                    L.bracket(
                        R, L.from_vector(phi * Y.to_vector())
                    ).to_vector()
                )
            )
            RHS = (
                2*alpha*(
                    L.bracket(L.bracket(R, X), Y).to_vector()
                    +
                    L.bracket(
                        L.bracket(R, L.from_vector(
                            phi * X.to_vector()
                            )
                        ),
                        L.from_vector(phi * Y.to_vector())
                    ).to_vector()
                )
            )[0]
            if debug:
                print("LHS =", LHS)
                print("RHS =", RHS)
            assert LHS == RHS, "g is not critical"
    print("g is critical!")

def Jordan_form_5d():
    """The (potentially) simplest example with ad_R
    _not_ diagonalizable: ad_R has two generalized
    eigenspaces, each of dimension two on which it acts
    by the matrix
    [ ±μ  ±ε ]
    [ 0   ±μ ]
    """
    alpha = matrix([1, 0, 0, 0, 0])
    g = identity_matrix(5)
    phi = matrix([
        [ 0, 0, 0, 0, 0],
        [ 0, 0, 0,-1, 0],
        [ 0, 0, 0, 0,-1],
        [ 0, 1, 0, 0, 0],
        [ 0, 0, 1, 0, 0]
        ])
    beta = g * phi

    check_compatible_structure(alpha, beta, phi, g)

    var('mu,epsilon')
    L = LieAlgebra(SR,
                    'R,v1,v2,w1,w2',
                    {
                        ('R', 'v1') : {'v1' : mu},
                        ('R', 'v2') : {'v1' : epsilon, 'v2' : mu},
    #                    ('R', 'v2') : {'v2' : mu},
                        ('R', 'w1') : {'w1' : -mu},
                        ('R', 'w2') : {'w1' : -epsilon, 'w2' : -mu},
    #                    ('R', 'w2') : {'w2' : -mu},
                    },
                   )
    R, v1, v2, w1, w2 = L.gens()

    check_alpha_R_invariant(L, alpha, beta, phi, g, R)
    check_beta_R_invariant(L, alpha, beta, phi, g, R)
    check_critical(L, alpha, beta, phi, g, R)
