#!/usr/bin/env julia

@enum sketchType subg=1 ros=2

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
	throw(ErrorException("Unimplemented operation!"))
end
