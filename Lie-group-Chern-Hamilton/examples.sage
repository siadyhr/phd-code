def non_ideal_zero_space():
    generators = "R,v1p,v1m,v2p,v2m"
    structure_coefficients = {
            ('R', 'v2p') : {'v2p' : 1},
            ('R', 'v2m') : {'v2p' : -1},
            ('v2p', 'v1p') : {'v2p' : 1},
            ('v2p', 'v1m') : {'v2p' : 1},
            ('v2m', 'v1p') : {'v2m' : 1},
            ('v2m', 'v1m') : {'v2m' : 1},
    }
    lie_algebra = LieAlgebra(QQ, generators, structure_coefficients)
    TestSuite(lie_algebra).run()

    R, v1p, v1m, v2p, v2m = lie_algebra.basis()
    h = lie_algebra.subalgebra([v1p,v1m])
    print("The subalgebra h")
    print(h)
    print("Is h an ideal of g?")
    print(h.is_ideal(lie_algebra))
    print(lie_algebra.representation('matrix'))

def non_ideal_zero_space_matrices():
    """Idea: The matrices representing the adjoint representation
    for the basis (R, v1pm, v2pm)
    v1pm have eigenvalue 0 for [R, -] while v2pm have eigenvalue pm 1.
    Also [v2pm, v1.] = v2pm (i.e. the v2. are -1-eigenvectors of the v1.)
    """
    R = diagonal_matrix([0,0,0,1,-1])
    v1p = matrix([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0,-1, 0],
        [0, 0, 0, 0,-1],
        ])
    v1m = v1p
    v2p = matrix([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [-1,1, 1, 0, 0],
        [0, 0, 0, 0, 0],
    ])
    v2m = matrix([
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0],
        [1, 1, 1, 0, 0],
        ])
    var('a b c d e')
    print(exp(a*R + b*v1p + c*v1m + d*v2p + e*v2m).simplify_full())

def double_eigenvalue():
    generators = "R,v1p,v1m,v2p,v2m"
    l = 1
    bp = 1
    bm = 1

    structure_coefficients = {
            ('R', 'v2p') : {'v2p' : 2*l},
            ('R', 'v2m') : {'v2m' : -2*l},
            ('R', 'v1p') : {'v1p' : l},
            ('R', 'v1m') : {'v1m' : -l},

            ('v2p', 'v2m') : {'R' : bp*bm/l},
            ('v2p', 'v1m') : {'v1p' : bp},
            ('v2m', 'v1p') : {'v1m' : bm},
    }
    lie_algebra = LieAlgebra(QQ, generators, structure_coefficients)
    TestSuite(lie_algebra).run()

double_eigenvalue()
