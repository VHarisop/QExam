module Sketch

	using Hadamard
	using StatsBase

	@enum sketchType subg=1 ros=2

	# low rank matrix generation
	lowrank(n, d, rank) = randn(n, rank) * randn(rank, d)

	"""
		fwht_no_scale!(y)

	Performs the fast Walsh-Hadamard transform, without scaling the resulting
	vector by ``1 / \sqrt{n}``. The vector is modified in place.
	"""
	function fwht_no_scale!(y)
		n = length(y)
		o = trailing_zeros(n) - 1
		m = i = n >> 1
		while i > 0
			for j = 0:m-1
				s = j + j >> o << o + 1
				t = s + i
				p, q = y[s], y[t]
				y[s], y[t] = p + q, p - q
			end
			i >>= 1; o -= 1
		end
	end


	"""
		fwht!(y, x)

	Performs the fast Walsh-Hadamard transform on `x`, scaling the resulting
	vector appropriately. The result is stored on `y`.
	"""
	function fwht!(y, x)
		n = length(x)
		ispow2(n) || throw(ArgumentError("Non power-of-2 HT length"))
		copy!(y, x)
		fwht_no_scale!(y)
		scale!(y, 1.0 / sqrt(n))
	end


	"""
		A_safe!(r, x, d)

	Computes the mapping ``x \\mapsto A x`` with ``A`` being a partial Hadamard
	matrix formed by keeping the first `length(x)` columns of ``H_m``. The
	result is stored in `r`. `d` plays the role of the diagonal matrix ``D``.
	"""
	function A_safe!(r, x, d)
		# compute D * x = d .* x
		fwht!(r, vcat(d .* x, zeros(length(r) - length(x))))
	end


	"""
		sketchmat_subg(width::Number, delta::Float64, n::Int, c0=3)

	Creates a sketching matrix for the subgaussian case, which requires ``m``
	to satisfy ``m \\geq \\left( \\frac{c_0 \\mathbb{W}(A \\mathcal{K})}{\\delta}
	\\right)^2``. ``c_0``, which can be optionally provided, is the constant
	appearing in the statement of the main Theorem.
	"""
	function sketchmat_subg(width::Number, delta::Float64, n::Int, c0=3)
		m = ceil(Int, c0 * (width / delta)^2)
		S = randn(m, n); return S
	end


	"""
		sketchmat_ros(width_s, width_r, delta, n, c0=3)

	Creates an operator that corresponds to multiplication with a sketching
	matrix coming from a randomized orthonormal system. `width_s` is the
	S-gaussian width and `width_r` is the Rademacher width of the tangent cone
	of the underlying problem. ``c_0``, which can be optionally provided, is the
	constant appearing in the statement of the main Theorem.
	"""
	function sketchmat_ros(width_s, width_r, delta, n, c0=3))
		lowstep = (c0 * (width_s^2 + log(n)) * (width_s^2)) / (delta^2)
		m = 2^(ceil(Int, log(lowstep)))  # Hadamard matrix must be power of 2
		D = rand([-1.0, 1.0], n)  # random diagonal vector
		Aop = (r, x) -> begin
			temp = zeros(m);
			A_safe!(temp, x, D);
			copy!(r, sample(temp, m, replace=true))
		end
		return Aop
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
			return sketchls_ros(A, x, delta)
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
		return xhat, xopt
	end


	function sketchls_ros(A, x, delta)
		throw(ErrorException("Unimplemented operation!"))
	end

end
