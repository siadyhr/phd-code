var('x1 x2 y1 y2')
alpha_sphere(x1, y1, x2, y2) = matrix([[-y1, x1, -y2, x2]])
dif_alpha_sphere(x1, y1, x2, y2) = matrix([
    [ 0, 2, 0, 0],
    [-2, 0, 0, 0],
    [ 0, 0, 0, 2],
    [ 0, 0,-2, 0],
])
Reeb_sphere(x1, y1, x2, y2) = vector([-y1, x1, -y2, x2])

print("Kerne af d alpha")
def sphere_tangent_projection(x1, y1, x2, y2):
    point = matrix([[x1, y1, x2, y2]]).T
    return identity_matrix(4) - point * point.T
print((
    sphere_tangent_projection(x1, y1, x2, y2).T
    *
    dif_alpha_sphere(x1, y1, x2, y2)
    *
    sphere_tangent_projection(x1, y1, x2, y2)
    ))

stereographic_projection(x1, y1, x2, y2) = (x1/(1-y2), y1/(1-y2), x2/(1-y2))
stereographic_projection_inverse(X1, Y1, X2) = vector(
        (
            2*X1,
            2*Y1,
            2*X2,
            (X1**2 + Y1**2 + X2**2 - 1)
        )
    )/(X1**2 + Y1**2 + X2**2 + 1)

alpha_projected(X1, Y1, X2) = (
        alpha_sphere(*stereographic_projection_inverse(X1, Y1, X2))
        *
        jacobian(
            stereographic_projection_inverse
            ,
            (X1, Y1, X2)
        )(X1, Y1, X2)
    )
print((1 + X1**2 + X2**2 + Y1**2)**2*alpha_projected(X1, Y1, X2).expand().simplify_full().factor())

Reeb_projected(X1, X2, Y2) = (
        jacobian(
            stereographic_projection,
            (x1, y1, x2, y2)
        )(
            *stereographic_projection_inverse(X1, X2, Y2)
        )
        *
        Reeb_sphere(*stereographic_projection_inverse(X1, X2, Y2))
    )

metric_projected(X1, Y1, X2) = (
        jacobian(
            stereographic_projection_inverse,
            (X1, Y1, X2)
        )(X1, Y1, X2).T
        *
        jacobian(
            stereographic_projection_inverse,
            (X1, Y1, X2)
        )(X1, Y1, X2)
    )

def exterior_derivative(one_form, variables):
    print(variables)
    print(12*"=")
    def out(variables):
        return variables[0]
    return out
    for i, entry in enumerate(one_form[0]):
        for j, variable in enumerate(variables):
            print(i, j, entry)
            print(variable)
            print(derivative(entry, variable))
            print(derivative(one_form[0, j], variables[i]))
            print("------")
            out(variables)[i, j] = derivative(entry, variable) - derivative(one_form[0, j], variables[i])
    return out

#dif_alpha_projected = exterior_derivative(alpha_projected, (X1, Y1, X2))
dif_alpha_projected(X1, Y1, X2) = (
        jacobian(stereographic_projection_inverse, (X1, Y1, X2))(X1, Y1, X2).T
        *
        dif_alpha_sphere(0,0,0,0)
        *
        jacobian(stereographic_projection_inverse, (X1, Y1, X2))(X1, Y1, X2)
        )
print((X1**2 + X2**2 + Y1**2 + 1)**3*dif_alpha_projected(X1, X2, Y1).simplify_full().factor())
print("Basis for dalpha efter projektion")
print(dif_alpha_projected(X1, Y1, X2).right_kernel().basis()[0].simplify_full())
R_candidate(X1, X2, Y1) = dif_alpha_projected(X1, Y1, X2).right_kernel().basis()[0]
print("Værdi af alpha")
print((alpha_projected(X1, Y1, X2) * R_candidate(X1, Y1, X2)).simplify_full())
print("Projektion af Reebfelt")
print(Reeb_projected(X1, Y1, X2).simplify_full())
print("Projektion af standardmetrikken")
print(metric_projected(X1, Y1, X2).simplify_full().factor())

#test_form(x, y) = matrix([[y, -x]])
#print(exterior_derivative(test_form, (x, y)))
