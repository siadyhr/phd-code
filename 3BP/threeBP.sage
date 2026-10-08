def V_constructor(mu):
    def V(x, y):
        return (
                (1-mu)/sqrt((x+mu)**2 + y**2)
                +
                mu/sqrt((x+mu-1)**2 + y**2)
            )
    return V

def Veff_constructor(mu):
    V = V_constructor(mu)
    Veff = (x**2 + y**2)/2 + V(x, y)
    return Veff

var('mu c x y q0')

Veff = Veff_constructor(mu)

print(
        (
            diff(Veff, x)(x=q0, y=0) #*(x-q0) + diff(Veff, y)*y
        )
    )
