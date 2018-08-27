#!/usr/bin/env julia

@enum sketchType subg=1 ros=2

# low rank matrix generation
lowrank(n, d, rank) = randn(n, rank) * randn(rank, d)


"""
	sketchmat_subg(width::Number, delta::Float64, n::Int)

Creates a sketching matrix for the subgaussian case, which requires ``m``
to satisfy ``m \\geq \\left( \\frac{c_0 \\mathbb{W}(A \\mathcal{K})}{\\delta}
\\right)^2``.
"""
function sketchmat_subg(width::Number, delta::Float64, n::Int)
	m = 2 * (width / delta)^2
	S = randn(m, n); return S
end


"""
"""
function sketchmat_ros(width_s::Number, width_r::Number, delta::Float64, n::Int)
	m = (2 * (width_s^2 + log(n)) * (width_s^2)) / (delta^2)
	throw(ErrorException("Unimplemented operation!"))
end


"""
	sketchls(A, x, stype::sketchType, delta)

Sketch and solve a least squares problem aiming for a ``\\delta`` approximation
under an appropriate sketch type.
"""
function sketchls(A, x, stype::sketchType, delta)
	if stype == subg
		return sketchls_subg(A, x, delta)
	elseif stype == ros
		throw(ErrorException("Unimplemented operation!"))
	end
end


"""
	sketchls_subg(A, x, delta)

Sketches a least squares problem without constraints using a subgaussian matrix
and solves the problem analytically. Returns the ratio
``\\frac{\\hat{x}}{x^{\\star}}``, where ``\\hat{x}`` is the solution to the
sketch problem and ``x^{\\star}`` is the solution to the original problem.
"""
function sketchls_subg(A, x, delta)
	n, d = size(A); r = rank(A); S = sketchmat_subg(sqrt(r), delta, n)
	Q = S * A; y = A * x; b = S * y
	xopt = A \ y; xhat = Q \ b
	return xhat / xopt
end
