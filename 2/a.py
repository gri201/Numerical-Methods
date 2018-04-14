import matplotlib.pyplot as plt
import numpy as np

def F(x):
	return np.cos(x)

def K(x):
	if x < 3. or x > 7.: return 1.
	elif x > 3. and x < 7.: return 3.

def Q(x):
	if x < 3. or x > 7.: return 3.
	elif x > 3. and x < 7.: return 1.

def f(n, N):
	return (np.sin((10. / N) * (n + 1./2)) - np.sin((10. / N) * (n - 1./2))) * N / 10. 

def b(n, N):
	return f(n - 1, N)

def a(n, N): 
	return K((10. / N) * (n + 1./2))

def q(n, N): 
	return (Q((10. / N) * (n - 1./2)) + Q((10. / N) * (n + 1./2))) / 2.

def alpha(n, N):
	return a(n - 1, N) / np.power(10. / N, 2) 

def beta(n, N):
	if n == 0: 
		return (K(0) * N / 10) + Q(0) * (10. / N) / 2
	else: 
		return (a(n - 1, N) + a(n, N)) / np.power(10. / N, 2) + q(n, N)

def gamma(n, N):
	if n == 0: 
		return K(0) / (10. / N)
	else: 
		return a(n, N) / np.power(10. / N, 2)

def Plot(c):
	plt.plot(c[:,0],c[:,1])
	plt.xlabel('x')
	plt.ylabel('y')
	plt.show()

def Solve(N):
	delta_1 = (10. / N) * F(0) / 2.
	phi_0 = -delta_1 / beta(0, N) 
	C_0 = gamma(0, N) / beta(0, N)

	phi = [0.0]
	C = [0.0]
	x = []
	y = []

	for i in range(0, N+1):
		x.append(i * 10. / N)

	phi.append(phi_0)
	C.append(C_0)

	for n in range(1, N-1):
		phi.append((alpha(n+1, N) * phi[-1] - b(n+1,N)) / (beta(n+1, N) - alpha(n+1, N) * C[-1]))
		C.append(gamma(n+1, N) / (beta(n+1, N) - alpha(n+1, N) * C[-1]))

	phi.append((alpha(N-1, N) * phi[-1] - b(N-1,N)) / (beta(N-1, N) - alpha(N-1, N) * C[-1]))

	y.append(0.)
	y.append(phi[N-1])
	for n in range(N-2, -1, -1):
		y.append(C[n] * y[-1] + phi[n])

	c = np.column_stack((np.array(x),np.array(list(reversed(y)))))

	h = 0.
	for element in c:
		print element[0], element[1], element[1] - h
		h = element[1]
	return c

def Error(N, c_1, c_2):
	err1 = np.fabs(c_1[0, 1] - c_2[0,1])
	err2 = np.fabs(c_1[N, 1] - c_2[2*N,1])
	err3 = 0

	for i in range(1, N):
		err3 += np.fabs(c_1[i, 1] - c_2[2 * i, 1])
	for i in range(1, N + 1):
		err3 += np.fabs((c_1[i,1 ] + c_1[i - 1, 1]) / 2 - c_2[2 * i - 1, 1])
	err3 = err3 * 10 / (2 * N)

	return err1 + err2 + err3


def RungeSolve(N, Eps):
	N_1 = N
	N_2 = 2 * N

	c_1 = Solve(N_1)
	c_2 = Solve(N_2)

	while Error(N_1, c_1, c_2) > Eps:
		N_1 = N_2
		N_2 *= 2

		c_1 = np.copy(c_2)
		c_2 = Solve(N_2)

	print 'Runge N: ', N_2
	return c_2

def microScope(c):
	print 'Microscope mode'
	a = input("Please, type a: ")
	b = input("Please, type b: ")

	d = np.column_stack(([], []))

	for element in c:
		if element[0] < b and element[0] > a:
			d = np.vstack((d, element))

	print 'Printing all nodes in interval [',a, ' , ',b,'] and it`s difference:'
	#print np.around(d, decimals=5)
	h = 0.
	for element in d:
		print element[0], element[1], element[1] - h
		h = element[1]

	Plot(d)



N = 1000
Eps = 0.01

#c = RungeSolve(N, Eps)
#Plot(c)

#N = input("Please, type custom N: ")
c = Solve(N)
Plot(c)

microScope(c)
#microScope(c)



