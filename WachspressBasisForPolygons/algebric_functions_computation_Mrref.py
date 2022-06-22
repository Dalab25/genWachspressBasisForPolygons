# -*- coding: utf-8 -*-
import sympy
import numpy as np
from algebric_functions import almostEqualRelativeAndAbs
from save_files_functions import save_M_rref_txt, save_pos_free_variables
from convert_data_into_variables import convert_data_into_free_pos_variables, convert_data_into_matrix_M_rref

def modify_M_to_Atild(order, pitch, bool_constraints=False, nSides=6):
	"""
	Function to change the matrix M_rref into a matrix A_tild and a vector b_tild
	"""
	unknowns_i = (order*(order + 1)/2) + (order -1)*(order*(order-1))/2
	nb_cols = int(nSides*unknowns_i)

	M_rref = convert_data_into_matrix_M_rref(order, pitch, bool_constraints, nSides)
	rank = np.linalg.matrix_rank(M_rref)
	#rank = M_rref.shape[1] - (order - 1)*4 - 1 #a  voir si on peut generaliser Ã  l'ordre k, 4 = nb de cotes "libres"
	A_tild = M_rref[:rank, :nb_cols]
	b_tild = M_rref[:rank,len(M_rref[1]) - 1]
	return A_tild, b_tild

def create_A_d_A_l(A_tild, order, pitch, nSides=6):
	"""
	Function to split the matrix A_tild into a matrix A_d with the dependent variables and A_l with the independent variables
	"""
	unknowns_i = (order*(order + 1)/2) + (order -1)*(order*(order-1))/2
	nb_cols = int(nSides*unknowns_i)
	pos_free_variables = convert_data_into_free_pos_variables(order, pitch, nSides)

	#matrix for the free variables
	A_l = []
	for i in range(len(pos_free_variables)):
		A_l.append(A_tild[:,pos_free_variables[i]])

	#matrix for the dependent variables
	A_d = []
	i = 0
	k = 0
	while (i != nb_cols - 1):
		j = pos_free_variables[k]
		if (i != j):
			A_d.append(A_tild[:,i])
			i+=1
		else:
			k+=1
			i+=1

	return A_d, A_l

def create_sol_x_l(order, pitch, bool_constraints=False, nSides=6):
	"""
	Function used to ensure that the linear systems set with and without the geometrical constraints are equivalent
	This function should be modified in order to change the independent variables
	"""

	pos_free_variables = convert_data_into_free_pos_variables(order, pitch, nSides)
	if(bool_constraints == True):
		solution = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(pitch) + '.txt')
	else:
		solution = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(pitch) + '.txt')

	x_l = []
	for i in pos_free_variables:
		x_l.append(solution[i])

	return x_l

def free_variables_equal_0(order, pitch, nSides=6):
	"""
	Fixing the free/independent variables to 0
	"""
	pos_free_variables = convert_data_into_free_pos_variables(order, pitch, nSides)
	x_l = np.zeros(np.shape(pos_free_variables))

	return x_l

def compute_x_d(order, pitch, bool_constraints, nSides=6):
	"""
	Computing x_d, for the moment, it is only done with x_l = [0,...,0]
	This function will be modified
	"""
	A_tild, b_tild = modify_M_to_Atild(order, pitch, bool_constraints, nSides)
	A_d, A_l = create_A_d_A_l(A_tild, order, pitch, bool_constraints, nSides)
	x_l = free_variables_equal_0(order, pitch) #WILL BE CHANGED
	x_d = np.linalg.solve(A_d, b_tild - np.dot(np.transpose(A_l),x_l))

	return x_d

def recreate_sol(x_l,x_d, order, pitch, bool_constraints, bool_free_variables_equal_0=False, nSides=6):
	"""
	Function to write the coefficients in the right order from x_l and x_d
	"""
	unknowns_i = (order*(order + 1)/2) + (order -1)*(order*(order-1))/2
	nb_unknowns = int(nSides*unknowns_i)

	pos_free_variables = convert_data_into_free_pos_variables(order, pitch, bool_constraints, nSides)
	x_l = np.array(x_l)

	sol = np.zeros(nb_unknowns)
	#In the aim of recreating the solution in the good order
	j = 0
	for i in range(nb_unknowns):
		k = 0
		for h in pos_free_variables:
			if(i == h):
				if(bool_free_variables_equal_0 == True):
					sol[i] = 0
				else:
					sol[i] = x_l[k]
				break
			else:
				k += 1
		if(k == len(pos_free_variables)):
			sol[i] = x_d[j]
			j += 1

	return sol

def recreate_sol_x_d_with_x_l_c_or_x_l_sc(order, pitch, bool_constraints, compare_solution=False, nSides=6):
	"""
	Computing the solutions this way
	x_d_c = A_d_c^{-1}[\tilde{b_c} - A_l_c x_l_wc]
	x_d_wc = A_d_wc^{-1}[\tilde{b_wc} - A_l_wc x_l_c]
	
	With:
	x_d_c (resp x_d_wc) = the dependent variables with (resp without) the constraints

	This function is useful to numerically ensure that the linear systems 
	built with an without the geometrical constraints are equivalent
	"""
	if(bool_constraints == True):
		sol = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(pitch) + '.txt')
	else:
		sol = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(pitch) + '.txt')

	A_tild, b_tild = modify_M_to_Atild(order, pitch, bool_constraints, nSides)
	A_d, A_l = create_A_d_A_l(A_tild, order, pitch, bool_constraints, nSides)

	coeffs_diff = []
	i = 0
	err = []

	''' algorithm to compare only x_d '''

	if(compare_solution == False):
		if bool_constraints == False:
			x_l_c = create_sol_x_l(order, pitch, not bool_constraints)
			x_d_sc = np.linalg.solve(A_d, b_tild - np.dot(np.transpose(A_l),x_l_c))

			x_d_sol = create_x_d_from_solution(order, pitch, not bool_constraints, nSides)

			#err = |x_d_sc - \tilde{x}_d_c|
			while (i!= len(x_d_sol)):
				err.append(x_d_sc[i] - x_d_sol[i])
				if(almostEqualRelativeAndAbs(x_d_sc[i], x_d_sol[i]) == False):
					coeffs_diff.append(i)
				i+=1

		else:
			x_l_sc = create_sol_x_l(order, not bool_constraints, nSides)
			x_d_c = np.linalg.solve(A_d, b_tild - np.dot(np.transpose(A_l),x_l_sc))
			x_d_sol = create_x_d_from_solution(order, not bool_constraints, nSides)

			#err = |x_d_c - \tilde{x}_d_sc|
			while (i!= len(x_d_sol)):
				err.append(x_d_c[i] - x_d_sol[i])
				if(almostEqualRelativeAndAbs(x_d_c[i], x_d_sol[i]) == False):
					coeffs_diff.append(i)
				i+=1

	''' algorithm to compare the whole solution '''
	#ATTENTION, values of the free variables has an impact
	if(compare_solution == True):
		x_l = create_sol_x_l(order, pitch, not bool_constraints, nSides)
		x_d = np.linalg.solve(A_d, b_tild - np.dot(np.transpose(A_l),x_l))
		solution = recreate_sol(x_l,x_d, order, bool_constraints, nSides)
		err = []
		while (i!= len(solution)):
			err.append(solution[i] - sol[i])
			if(almostEqualRelativeAndAbs(solution[i], sol[i]) == False):
				coeffs_diff.append(i)
			i+=1

	if(len(coeffs_diff) == 0):
		print("Identical solutions")

	print(np.linalg.norm(err))

def create_x_d_from_solution(order,pitch, bool_constraints, nSides=6):
	"""
	Function to get the dependent coefficients from the least squares solution
	"""
	if(bool_constraints == True):
		sol = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(pitch) + '.txt')
	else:
		sol = np.loadtxt('./results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(pitch) + '.txt')

	pos_free_variables = convert_data_into_free_pos_variables(order, pitch, bool_constraints, nSides)

	x_d = []
	for i in range(len(sol)):
		if(pos_free_variables != []):
			for k in pos_free_variables:
				if(i == k):
					break
			if(i != k):
				x_d.append(sol[i])
		else:
			x_d.append(sol[i])

	return x_d

def computationSolWithFreeVariables0(order, p, boolConstraints, nSides=6):
	"""
	Computing the coefficients with the free variables fixed to 0
	"""
	if(boolConstraints):
		sol = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_constraints' + str(p) + '.txt')
	else:
		sol = np.loadtxt('../results/polygon_' + str(nSides) + 'sides/order' + str(order) + '/symbolic_coefficients/solution_order_' + str(order) + '_no_constraints' + str(p) + '.txt')

	freePosVariables = convert_data_into_free_pos_variables(order, p, nSides)
	A_tild, b_tild = modify_M_to_Atild(order, p, boolConstraints, nSides)
	A_d, _ = create_A_d_A_l(A_tild, order, p, nSides)
	x_d = np.linalg.solve(A_d, b_tild) #- np.dot(np.transpose(A_l),x_l))

	solution = np.zeros(len(sol))
	for i,freePosVariableI in enumerate(freePosVariables):
		solution[freePosVariableI] = 0

	counter = 0
	for i in range(len(solution)):
		if(i not in freePosVariables):
			solution[i] = x_d[counter]
			counter += 1
	np.savetxt("../results/polygon_" + str(nSides) + "sides/order" + str(order) + "/symbolic_coefficients/solution_order_" + str(order) + "_pitch_" + str(p) + "_free_variables_0.txt", solution)

def rref(B, tol=1e-10):
	"""
	Function rref taken from:
	https://gist.github.com/sgsfak/77a1c08ac8a9b0af77393b24e44c9547
	"""
	A = B.copy()
	rows, cols = A.shape
	r = 0
	pivots_pos = []
	row_exchanges = np.arange(rows)
	for c in range(cols):
		## Find the pivot row:
		pivot = np.argmax (np.abs (A[r:rows,c])) + r
		m = np.abs(A[pivot, c])
		if m <= tol:
			## Skip column c, making sure the approximately zero terms are
			## actually zero.
			A[r:rows, c] = np.zeros(rows-r)
		else:
			## keep track of bound variables
			pivots_pos.append((r,c))
			if pivot != r:
				## Swap current row and pivot row
				A[[pivot, r], c:cols] = A[[r, pivot], c:cols]
				row_exchanges[[pivot,r]] = row_exchanges[[r,pivot]]
			## Normalize pivot row
			A[r, c:cols] = A[r, c:cols] / A[r, c];
			## Eliminate the current column
			v = A[r, c:cols]
			## Above (before row r):
			if r > 0:
				ridx_above = np.arange(r)
				A[ridx_above, c:cols] = A[ridx_above, c:cols] - np.outer(v, A[ridx_above, c]).T
			## Below (after row r):
			if r < rows-1:
				ridx_below = np.arange(r+1,rows)
				A[ridx_below, c:cols] = A[ridx_below, c:cols] - np.outer(v, A[ridx_below, c]).T
			r += 1
		## Check if done
		if r == rows:
			break;
	return (A, pivots_pos, row_exchanges)


def compute_M_rref(A, B, constraints, order, pitch, save_M_rref, save_pos_free_variables, method_rref=0, nSides=6):
	#converting B in a good format
	b = np.zeros((len(B),1))
	for i in range(len(B)):
		b[i] = np.array([B[i]])

	M = np.hstack([A,b])
	if(method_rref != 0):
		print("Sympy rref method")
		M = sympy.Matrix(M)
		M_rref,nb_pivots = M.rref()
		print("Vector containg the pivot lines", nb_pivots)
		print("Rank M_rref = ",len(nb_pivots))
	else:
		print("Local rref method")
		M_rref,piv,_ = rref(M)
		pos_piv = get_pos_piv(piv)
		get_pos_free_variables(M,pos_piv,order, pitch, save_pos_free_variables, nSides)

	if(save_M_rref == True):
		save_M_rref_txt(order, pitch, M_rref, constraints, nSides)

	return M_rref

def get_pos_piv(piv):
	"""
	Get the placement of the pivots
	"""
	pos_piv = []
	for i in range(len(piv)):
		pos_piv.append(piv[i][1])
	return pos_piv

def get_pos_free_variables(M,pos_piv,order, pitch, save_free_variables, nSides):
	"""
	Function to get the placement of the free variables in the rref matrix
	"""
	pos_free_variables = []
	v = []
	for k in range(M.shape[1]):
			v.append(k)
	i = 0
	j = 0
	while (i != M.shape[1]):
		if(v[i] != pos_piv[j]):
			diff = pos_piv[j] - v[i]
			if (diff >= 0):
				for k in range(diff):
					pos_free_variables.append(v[i])
					i+=1
			else:
				pos_free_variables.append(v[i])
				i+=1
		else:
			i+=1
			if(j != len(pos_piv) - 1):
				j+=1
	if(pos_free_variables != []):
		if (pos_free_variables[len(pos_free_variables) - 1] == M.shape[1] - 1):
			pos_free_variables.pop()

	if(save_free_variables == True):
		print(pos_free_variables)
		save_pos_free_variables(order, pitch, pos_free_variables, nSides)
	return pos_free_variables
