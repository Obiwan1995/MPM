#!/usr/bin/python
# -*- coding: utf-8 -*-

from math import *
import numpy as np
import matplotlib.pyplot as pyplot

def f(A, xt, B, ut):
	return A.dot(xt) + B.dot(ut)

def g(A, pT):
	return np.transpose((-1)*A).dot(pT)

def changement_coord_x1(x1, x2, n, delta_t):
	return x1*cos(omega*n*delta_t) - x2*sin(omega*n*delta_t) + cos(omega*n*delta_t)

def changement_coord_x2(x1, x2, n, delta_t):
	return x1*sin(omega*n*delta_t) + x2*cos(omega*n*delta_t) + sin(omega*n*delta_t)

def resolution_systeme_fusee(A, B, x0, ut, delta_t):
	x = x0
	i = 0
	while i < n:
		val = f(A, x, B, ut)
		x = x + delta_t * val
		tab_x.append([x[0,0], x[2,0]])
		tab_x.append([changement_coord_x1(x[0,0], x[2,0], i, delta_t), x[1,0], changement_coord_x2(x[0,0], x[2,0], i, delta_t), x[3,0]])
		i += 1
	return x

def resolution_systeme_fusee2(A, B, x0, ut, delta_t):
	x=x0
	x2=0
	i=0
	while i < n:
		val = f(A, x, B, ut)
		val2 = f(A, x2, B, ut)
		x2=x+(delta_t/2)*val
		x=x+delta_t*val2
		#tab_x.append([x[0,0], x[2,0]])
		tab_x.append([changement_coord_x1(x[0,0], x[2,0], i, delta_t), x[1,0], changement_coord_x2(x[0,0], x[2,0], i, delta_t), x[3,0]])
		i+=1
	return x

def resolution_systeme_retrograde(A, xT, delta_t):
	t = 0
	p = xT
	i = 0
	while i < n:
		val = g(A, p)
		p = p - delta_t * val
		#t = t + delta_t
		i += 1
	return p

def resolution_systeme_retrograde2(A, xT, delta_t):
	t = 0
	p = xT
	p2 = 0
	i = 0
	while i < n:
		val = g(A, p)
		val2 = g(A,p2)
		p2=p + (delta_t/2)*val
		p = p - delta_t * val2
		#t = t + delta_t
		i += 1
	return p

def gradient_J(epsilon, u, B, p):
	return epsilon*u+np.transpose(B).dot(p)

def methode_gradient(epsilon, B, p, u0, ro):
	u = u0
	t = 0
	i = 0
	while i < 100:
		u = u - ro * gradient_J(epsilon, u, B, p)
		t = t + delta_t
		i += 1
	return u

#main
tab_x = []
T = 1
omega = 2*pi / T
A = np.matrix([	[0,1,0,0], 
				[3*omega*omega,0,0, 2*omega],
				[0,0,0,1],
				[0,-2*omega,0,0,] ]);

B = np.matrix([	[0,0], 
				[1,0], 
				[0,0], 
				[0,1]]);
x = np.transpose(np.atleast_2d([1,0,0,0]))
u = np.transpose(np.atleast_2d([0,0]))
delta_t = 0.001
n = T/delta_t
epsilon = 0.00000001
ro = 0.01
for i in range(0,50):
	x = resolution_systeme_fusee2(A, B, x, u, delta_t)
	p = resolution_systeme_retrograde2(A, x, delta_t)
	u = methode_gradient(epsilon, B, p, u, ro)

x1 = []
x2 = []
y1 = []
y2 = []
for i in range(0, len(tab_x)):
	x1.append(changement_coord_x1(tab_x[i][0], tab_x[i][1], i, delta_t))
	x2.append(changement_coord_x2(tab_x[i][0], tab_x[i][1], i, delta_t))
	#x1.append(tab_x[i][0])
	#x2.append(tab_x[i][1])
	y1.append(10**27*cos(omega*i*delta_t))
	y2.append(10**27*sin(omega*i*delta_t))
	
pyplot.plot(x1, x2, "r")
pyplot.plot(y1, y2, "b")
pyplot.axis("equal")

pyplot.show()

# Euler explicite
# x(n+1) = x(n) + delta(t)f(x(n))

# Euleur amélioré
# x(n+1) = x(n) + delta(t)*f(x(n) + 1/2)
# x(n+1/2) = x(n) + delta(t)/2 * f(x(n))

# Méthode du gradient
# u(n+1) = u(n) - p * gradient(J) avec p = pas de la méthode