# 3corposC

Final Project of my Computational Physics Course solving the 3 Body Problem for Solar System-ish cenarios.

The 2dstruc.h is where i've defined the structure for a 2-dimensional vector, it also includes a norm(vec) function to
get the norm of the vector into that format, and the most important thing, a Runge kutta of fourth order for second ordem and bi dimensional EDO's (Specially usefull to solve Newton's Equation in 2 Dimensions)

In the 2corpos.c I solved the 2 Body problem to test if my Runge Kutta program was correct - since I already knew the answer
I expected. After debbuging the Runge kutta I followed to:

3crps.c which uses my program of the Runge Kutta to solve the 3 body problem for real, the way i'm returning the solution is pretty confunsing, did that wrong, but was on the dead line, it returns a 3 dimension array, first indice is the body, second indice is: 0 - x position, 1 - y position, 2 - x speed, 3 - y speed.

PS: Unities I optimized for: Distance: Astronomical units (AU), Speed: AU/year, time: years , massa: Solar MAsses.
