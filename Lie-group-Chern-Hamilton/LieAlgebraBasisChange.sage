"""This script calculates basis changes of Lie
algebras and the pushforward of a cosymplectic
and metric structure under this change.

NB: Given a Lie bracket on generators (e_i) and
a linear map L, it constructs a new Lie algebra
with generators f_i := L e_i. The Lie bracket
is [f_i, f_j] := [L e_i, L e_j] (considering
L as a map of the original Lie algebra).

If one wants to think of L as a linear map
L : g -> h, where
- g has
    - basis e_i
    - bracket [e_i, e_j]
- h has
    - basis f_i
    - bracket [f_i, f_j] := [L e_i, L e_j]
one should instead use L^{-1}

"""
var('l l1 l2 m2 m3 a b c bp bm')

structure_coefficients = {
        'abelian-lambda1=0' : {
               ('e5', 'e2') : {'e2' : l},
               ('e5', 'e3') : {'e3' : -l},
           },
        'abelian-lambda1!=0' : {
               ('e5', 'e1') : {'e1' : l1},
               ('e5', 'e2') : {'e2' : l2},
               ('e5', 'e3') : {'e3' : -l2},
               ('e5', 'e4') : {'e4' : -l1},
           },
        'sol+R' : {
               ('e5', 'e2') : {'e2' : l},
               ('e5', 'e3') : {'e3' : -l},

               ('e1', 'e2') : {'e2' : 1},
               ('e1', 'e3') : {'e3' : -1},
            },
        'h3+R' : {
               ('e1', 'e2') : {'e3' : 1},

               ('e5', 'e1') : {'e1' : 2*l},
               ('e5', 'e2') : {'e2' : -l},
               ('e5', 'e3') : {'e3' : l},
               ('e5', 'e4') : {'e4' : -2*l},
            },
        'h3+R-raw' : {
               ('e1', 'e2') : {'e3' : b},

               ('e5', 'e1') : {'e1' : 2*l},
               ('e5', 'e2') : {'e2' : -l},
               ('e5', 'e3') : {'e3' : l},
               ('e5', 'e4') : {'e4' : -2*l},
            },

        'non-cosymplectic-lambda1=0' : {
            # Basis (R2, v+, v-, R3, R)
            ('e5', 'e2') : {'e2' : l},
            ('e5', 'e3') : {'e3' : -l},

            ('e1', 'e2') : {'e2' : m2},
            ('e1', 'e3') : {'e3' : -m2},
            
            ('e4', 'e2') : {'e2' : m3},
            ('e4', 'e3') : {'e3' : -m3},

            ('e2', 'e3') : {'e5' : a},
            },
        'non-cosymplectic-lambda1!=0' : {
            # Basis (v2+, v1+, v1-, v2-, R)
               ('e1', 'e3') : {'e2' : bp},
               ('e4', 'e2') : {'e3' : bm},
               ('e1', 'e4') : {'e5' : a},

               ('e5', 'e1') : {'e1' : 2*l},
               ('e5', 'e2') : {'e2' : l},
               ('e5', 'e3') : {'e3' : -l},
               ('e5', 'e4') : {'e4' : -2*l},
            }
}

basis_changes = {
        'abelian-lambda1=0' : matrix([
            [0,     0,  0,  0,  1],
            [0,     1,  0,  0,  0],
            [0,     0,  1,  0,  0],
            [0,     0,  0,  1,  0],
            [1/l,   0,  0,  0,  0],
        ]),
        'abelian-lambda1!=0' : matrix([
            [0, 0,  1,  0,  0],
            [0, 0,  0,  1,  0],
            [-1,0,  0,  0,  0],
            [0,-1,  0,  0,  0],
            [0, 0,  0,  0,  1/l2],
        ]),
        'sol+R' : matrix([
            [1/2,       0,  0,  0,  1/2],
            [0,         1,  0,  0,  0],
            [0,         0,  1,  0,  0],
            [0,         0,  0,  1,  0],
            [1/(2*l),   0,  0,  0,  -1/(2*l)]
        ]),
        'h3+R' : matrix([
            [0, 0, 1, 0, 0],
            [0,-1, 0, 0, 0],
            [1, 0, 0, 0, 0],
            [0, 0, 0, 1, 0],
            [0, 0, 0, 0, 1/l]
        ]),
        'h3+R-raw' : matrix([
            [0, 0, 1, 0,    0],
            [0,-1, 0, 0,    0],
            [b, 0, 0, 0,    0],
            [0, 0, 0, 1/b,  0],
            [0, 0, 0, 0,    1/l]
        ]),
        'non-cosymplectic-lambda1=0' : matrix([
            [0,     0,              0,          l,      0],         # R2
            [0,     1/(l*a),        0,          0,      0],         # v+
            [0,     0,              1,          0,      0],         # v-
            [0,     0,              0,          0,      l],         # R3
            [1/l,   0,              0,          -m2,  -m3],     # R

            # Cols: (sl 2) + R^2
        ]),
        'non-cosymplectic-lambda1!=0' : matrix([
            [0,     1/bp,
                            0,      0,          0],
            [0,     0,      0,      1,  0],
            [0,     0,      0,      0,          1],
            [0,     0,      1/bm,
                                    0,          0],
            [1/l,
                    0,      0,      0,          0]
             # Rows: v2+, v1+, v1-, v2-, R
             # Cols: (sl 2) + R^2
        ]),
}

metrics = {
        'abelian-lambda1=0' : identity_matrix(5),
        'abelian-lambda1!=0' : identity_matrix(5),
        'sol+R' : diagonal_matrix([c**2, 1, 1, c**(-2), 1]),
        'h3+R' : diagonal_matrix([c**(-2), 1, 1, c**(2), 1]),
        'h3+R-raw' : identity_matrix(5),
        'non-cosymplectic-lambda1=0' : identity_matrix(5),
        'non-cosymplectic-lambda1!=0' : identity_matrix(5)
        }

def converter(model):
    print("Model:", model)
    lie_algebra = LieAlgebra(SR,
         'e1,e2,e3,e4,e5',
         structure_coefficients[model]
    )
    basis_change = basis_changes[model]
    print("Basis change:")
    print(basis_change)
    print("Basis change inverse:")
    print(1/basis_change)
    fs = [
            lie_algebra.from_vector(basis_change * ei.to_vector())
            for ei in lie_algebra.basis()
        ]
    print(fs)
    new_lie_algebra = LieAlgebra(SR, 'f1,f2,f3,f4,f5',
        {
            ('f%s' % i, 'f%s' % j) : dict(zip(
                'f1,f2,f3,f4,f5'.split(","),
                (1/basis_change) * (fi.bracket(fj)).to_vector()
                ))
            for j, fj in enumerate(fs, 1)
            for i, fi in enumerate(fs, 1)
        }
     )

    print(new_lie_algebra.structure_coefficients())

    metric_matrix = basis_change.T * metrics[model] * basis_change

    print("Metric in the new basis:")
    print(metric_matrix)

    print("Cosymplectic 1-form in the new basis:")
    cosymplectic_1_form_matrix = basis_change[4,:][0]
    print(" + ".join(
        "%s e^%s" % (a, i+1)
        for i, a in enumerate(cosymplectic_1_form_matrix)
        if a != 0
        ))

    print("Cosymplectic 2-form in the new basis:")
    beta_0 = matrix([
                [ 0, 0, 0, 1, 0], 
                [ 0, 0, 1, 0, 0], 
                [ 0,-1, 0, 0, 0], 
                [-1, 0, 0, 0, 0], 
                [ 0, 0, 0, 0, 0]
            ]) # e^14 + e^23
#    beta_0 = matrix([
#                [ 0, 1, 0, 0, 0],
#                [-1, 0, 0, 0, 0],
#                [ 0, 0, 0, 1, 0],
#                [ 0, 0,-1, 0, 0],
#                [ 0, 0, 0, 0, 0]
#            ]) # e^12 + e^34
    cosymplectic_2_form_matrix = (
            basis_change.T
            *
            beta_0
            *
            basis_change
        )
    #print(cosymplectic_2_form_matrix)
    print(" + ".join(
        [
            "%s e^%s%s" % (a_ij, i+1, j+1)
            for i, row in enumerate(cosymplectic_2_form_matrix)
            for j, a_ij in enumerate(row[i+1:], i+1) if a_ij != 0
        ]
    ))

models = [
        'abelian-lambda1=0',
        'abelian-lambda1!=0',
        'sol+R',
        'h3+R',
        'h3+R-raw',
        'non-cosymplectic-lambda1=0',
        'non-cosymplectic-lambda1!=0'
        ]
for L in basis_changes.items():
    break
    print(L[0])
    print(L[1])
    print()
    print(1/L[1])
    print()

converter(models[0])
