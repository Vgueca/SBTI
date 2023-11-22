from sympy import groebner
from sympy.abc import x,y,a,b
from obtener_base import Tjurina_generator

F = [y**5 - x**7 + 1**x**3*y**3 + 1*x**4*y**4]

print(groebner(F,x,y,domain="CC",order="grevlex"))