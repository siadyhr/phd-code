import functools
load("LieAlgebraBasisChange.sage")

def Riem_constructor(Lie_algebra, g):
    """Compute the Riemannian curvature tensor
    using formula for the Christoffel symbols.
    The implementation only works when the

    an orthogonal basis to work.

    Parameters
    ----------
    Lie_algebra : LieAlgebra with n generators
    g : n×n matrix
        Represents the inner product in the standard
        basis `Lie_algebra.basis()` of Lie_algebra

    Returns
    -------
    Riem : function
        (X, Y, Z, W) -> R(X, Y, Z, W)
        Computes the Riemannian curvature 4-tensor.
        The input should be elements in Lie_algebra.

    Internal helper functions:

    Christoffel : function
        (i, j, k) -> Gamma_ij^k
        Needs `Lie_algebra.basis()` to be orthogonal
        wrt `g` to work

    c : function
        (i, j, k) -> c_ij^k
        Extracts the structure coefficients
        of Lie_algebra:
        [ei, ej] = sum_k c_ij^k ek, where
        ei etc. are the elements of `Lie_algebra.basis()`

    Riem_coeff : function
        (i, j, k, l) -> R_ijkl
        Computes the components of the Riemannian
        4-tensor in the ORTHOGONAL basis
        `Lie_algebra.basis()`

    Helper functions are cached using `functools.cache`
    to speed up computations
    """

    # Ensure that Lie_algebra.basis() is orthogonal
    assert (
            metrics[models[model_id]]
            -
            diagonal_matrix(metrics[models[model_id]].diagonal())
            ==
            0
        )
    basis = list(L.basis())

    @functools.cache
    def Christoffel(i,j,k):
        """Calculate Christoffel symbol Gamma_ij^k, i.e.,
        nabla_(ei) ej = sum_k Gamma_ij^k ek.

        Uses formula from Milnor 1976,
        "Curvatures of left invariant metrics on lie groups":
        g(nabla_X Y, Z)
        =
        (1/2) (
            g([X, Y], Z)
            -
            g([Y, Z], X)
            +
            g([Z, X], Y)
        )
        When the basis is orthogonal, this implies
        Gamma_ij^k = (1/g(ek, ek)) * g(nabla_(ei) ej, ek)
        """

        # basis[i] are elements in Lie_algebra
        # .to_vector() lets us multiply them with
        # the matrix g
        return (1/2)*(1/g[k, k]) * (
                (basis[i].bracket(basis[j])).to_vector() * g * basis[k].to_vector()
                -
                (basis[j].bracket(basis[k])).to_vector() * g * basis[i].to_vector()
                +
                (basis[k].bracket(basis[i])).to_vector() * g * basis[j].to_vector()
            )

    @functools.cache
    def c(i, j, k):
        # Utility function to calculate the
        # structure coefficients of the Lie algebra
        # L.basis().keys()[k] gives the string representing
        # the k'th basis element
        return (basis[i].bracket(basis[j])).coefficient(L.basis().keys()[k])

    @functools.cache
    def Riem_coeff(i, j, k, l):
        """Calculate components of the
        Riemannian 4-tensor in the basis
        `Lie_algebra.basis()`
        """
        return sum(
            sum(
                g[l, m] * (
                    Christoffel(j, k, m) * Christoffel(i, m, n)
                    -
                    Christoffel(i, k, m) * Christoffel(j, m, n)
                    -
                    c(i, j, m) * Christoffel(m, k, n)
                )
                for m, _ in enumerate(L.basis())
            )
            for n, _ in enumerate(L.basis())
        )

    def Riem(X, Y, Z, W):
        """Calculate Riemannian 4-tensor on elements
        from `Lie_algebra`
        """
        return sum(
                sum(
                    sum(
                        sum(
                            xi * yj * zk * wl * Riem_coeff(i, j, k, l)
                            for i, xi in enumerate(X.to_vector()) if xi != 0
                        )
                        for j, yj in enumerate(Y.to_vector()) if yj != 0
                    )
                    for k, zk in enumerate(Z.to_vector()) if zk != 0
                )
                for l, wl in enumerate(W.to_vector()) if wl != 0
            )
        # The checks `if xi != 0` etc. avoids
        # unnecessary calls to Riem_coeff
    return Riem

def Ricci_constructor2(Lie_alg, g):
    """Compute the 2-Ricci tensor on
    elements from `Lie_alg`
    """
    Riem = Riem_constructor(Lie_alg, g)
    def Ricci2(X, Y):
        return sum(
                Riem(ek, X, Y, ek)/g[k,k]
                for k, ek in enumerate(Lie_alg.basis())
                )
    return Ricci2

def Ricci_constructor11(Lie_alg, g):
    """Compute the (1, 1)-Ricci tensor on
    elements from `Lie_alg`
    """
    Ric2 = Ricci_constructor2(Lie_alg, g)
    @functools.cache
    def Ricci11(X):
        return sum(
                (1/g[j,j]) * Ric2(X, ej) * ej
                for j, ej in enumerate(Lie_alg.basis())
                )
    return Ricci11

def soliton_equations(Lie_alg, g):
    """ (Lie_alg, g) is a Ricci-soliton
    <=>
    There exists a constant C
    st. T is a derivation
    """
    var('C')
    Ric11 = Ricci_constructor11(Lie_alg, g)
    def T(X):
        return Ric11(X) - C*X

    for ei in L.basis():
        for ej in L.basis():
            # eq vanishes for all ei, ej <=> T is a derivation
            eq = T(ei.bracket(ej)) - (T(ei)).bracket(ej) - ei.bracket(T(ej))
            if eq != 0:
                print("Check (%s, %s) non-trivial:" % (ei, ej))
                print(eq)
                print("Non-trivial component equations:")
                for subeq in eq.to_vector():
                    # eq is a vector equation; split
                    # into one equation for each
                    # component
                    if subeq == 0:
                        continue
                    print("\t", subeq)
                    # Too substitute variables with values,
                    # use .substitute, e.g.
#                    print("\t", subeq.substitute({m2 : 0}))
                    print("\t <=>", solve(subeq==0, C)[0])

# Model Lie algebras as defined in
# `LieAlgebraBasisChange.sage`
model_id = 0
L = LieAlgebra(SR,
               'e1,e2,e3,e4,e5',
               structure_coefficients[models[model_id]]
           )
soliton_equations(L, metrics[models[model_id]])
