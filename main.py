#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import *
import numpy as np
import matplotlib.pyplot as pyplot

def f(xt, ut):
	return np.dot(A, xt) + np.dot(B, ut)

def g(pT):
	return -np.dot(A.T, pT)

def changement_coord_x1(x1, x2, n):
	return x1*cos(omega*n*delta_t) - x2*sin(omega*n*delta_t) + cos(omega*n*delta_t)

def changement_coord_x2(x1, x2, n):
	return x1*sin(omega*n*delta_t) + x2*cos(omega*n*delta_t) + sin(omega*n*delta_t)

def resolution_systeme_fusee(x, u):
	tab_x = [x]
	for i in range(n):
		val = f(x, u)
		x = x + delta_t * val
		tab_x.append(x)
	return tab_x

def resolution_systeme_fusee2(x, u):
	tab_x = [x]
	for i in range(n):
		val = f(x, u)
		xn12 = x+(delta_t/2)*val
		val2 = f(xn12, u)
		x = x+delta_t*val2
		tab_x.append(x)
	return tab_x

def resolution_systeme_retrograde(p):
	tab_p = [p]
	for i in range(n):
		val = g(p)
		p = p - delta_t * val
		tab_p.append(p)
	return tab_p

def resolution_systeme_retrograde2(p):
	tab_p = [p]
	for i in range(n):
		val = g(p)
		p2=p - (delta_t/2)*val
		val2 = g(p2)
		p = p - delta_t * val2
		tab_p.append(p)
	return tab_p

def gradient_J(u, p):
	return epsilon*u+np.dot(B.T,p)

##### MAIN #####
T = 1.0
omega = (2.0*pi) / T
A = np.matrix([	[0.,1.0,0.,0.], 
				[3.*omega**2.,0.,0., 2.*omega],
				[0.,0.,0.,1.0],
				[0.,-2.*omega,0.,0.] ]);

B = np.matrix([	[0.,0.], 
				[1.0,0.], 
				[0.,0.], 
				[0.,1.0]]);
x0 = np.transpose(np.atleast_2d([0.1,0.,0.,0.]))
u0 = np.transpose(np.atleast_2d([0.,0.]))
n = 500
delta_t = float(T/n)
epsilon = 0.001
ro = 0.03

# CALCUL
for j in range(n):
	tab_u = [u0]
	tab_x = resolution_systeme_fusee(x0, u0)
	tab_p = resolution_systeme_retrograde(tab_x[n])

	u0 = u0 - ro * gradient_J(u0, tab_p[n-j])
	tab_u.append(u0)

x1 = []
x2 = []
y1 = []
y2 = []

for i in range(n):
	x1.append(changement_coord_x1(tab_x[i][0][0], tab_x[i][2][0], i))
	x2.append(changement_coord_x2(tab_x[i][0][0], tab_x[i][2][0], i))
	y1.append(cos(omega*i*delta_t))
	y2.append(sin(omega*i*delta_t))
	
pyplot.plot(x1, x2, "r")
pyplot.plot(y1, y2, "b")
#pyplot.axis("equal")

pyplot.show()

# Euler explicite
# x(n+1) = x(n) + delta(t)f(x(n))

# Euler amélioré
# x(n+1) = x(n) + delta(t)*f(x(n + 1/2))
# x(n+1/2) = x(n) + delta(t)/2 * f(x(n))

# Méthode du gradient
# u(n+1) = u(n) - p * gradient(J) avec p = pas de la méthode