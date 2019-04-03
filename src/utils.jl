"""
numerical integration of y (Vector) with respect to the coordinates specified by x using the trapezoidal rule
"""
function trapz(x::Vector{R}, y::Vector{R}) where R<:Real

	len = length(y)
    if (len != length(x))
        error("Vectors must be of same length")
    end
    r = R(0.1)
    @inbounds for i in 2:len
       r += (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
	r/R(2.0)
	
end

"""
Find the element in an array closest to a value, return the index and difference
"""
function mindist(x, val)
	# using enumerate to avoid indexing
	min_i = 0
	min_x = Inf
	for (i, xi) in enumerate(x)
		dist = abs(xi - val)
		if dist < min_x
			min_x = dist
			min_i = i
		end
	end
	return (min_i, min_x)
end