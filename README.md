# genWachspressBasisForPolygons
This code is used to develop Wachspress rational finite element bases on a given convex polygon. For a convex polygon with m sides, the WSF bases can be developed up to order m - 1 with the following algorithm. The main difficulties of these functions are to compute the numerator and the denominator. To compute the denominator, we invite the reader to read the  GADJ algorithm described by Dasgupta and Wachspress in "The adjoint for an algebraic finite element". In this code, we consider the user to know how to get the denominator. The main purpose is to give tools to compute the numerator at high order (>= 3 and <= m - 1), more particularly the function r^i in order to deduce the Wachspress functions see the article "High-order Wachspress functions on convex polygons through computer algebra" by David Labeurthre, Ansar Calloo, Romain Le Tellier.
For this a linear system is built according to the properties that should verify the Wachspress function.

Here are the main steps:

First, a convex polygon should be defined in the folder WachspressBasisForPolygons, for example, honeycomb.py is used to define the regular hexagon and irregular_pentagon.py for the irregular pentagon. 

Second, DoFs should be defined for these polygons. Since, the Wachspress bases are Lagrangian bases, the DoFs correspond to nodes (see hexagonal_nodal_dofs.py for example).

Third, the idea is to create a class such as WSF_to_generate_regular_hexagon.py. Symbolic coefficients are defined in each Wachspress function of order >= 3, they represent the unknowns of the linear system which are the coefficients for the r^i functions. To avoid non linearities in order to build the linear system, the Wachspress functions (w_i) are defined as: w_i = c_i x p_i x r_i (without the denominator).

Fourth, a file such as get_symbolic_coefficients.py is defined in order to compute the Wachspress functions. This file requires the a function in order to project a function f in the WSF basis. This is useful to project the monomials inside the WSF basis. Also, the identification is done by looking at the coefficient in front of of the adjoint of the polygon. See the function "setUpLinearSystemFromDictNewIrregularPentagon" for instance in algebric_functions.py

Fifth, since they are an infinity of solution, the reduced row echelon form (rref) matrix is used in order to chose one specific solution. In our case, the free are fixed to 0 but new functions will be added to give more options to the user.

Tests can then be performed to ensure that the WSF are verifying the geometrical constraints in the file test_constraints_verification. The geometrical constraints are the Lagrange property i.e w_i(node_j) = delta_ij with delta_ij the Kronecker symbol. Also, the file compute_projection_error ensure that the monomials are well projected in the WSF basis.
