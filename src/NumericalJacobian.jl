"""Differential-state names used when assembling a numerical Jacobian pattern."""
function jacobian_state_names end

jacobian_state_names(model) = throw(ArgumentError(
	"no Jacobian state layout is registered for model $(nameof(model))"))

function _compiled_state_indices(compiled, component, names)
	indices = Int[]
	for name in names
		hasproperty(component, name) || continue
		index = SII.variable_index(compiled, getproperty(component, name))
		index === nothing || push!(indices, index)
	end
	return unique(indices)
end

function _append_dense_block!(rows, cols, row_indices, column_indices)
	for column in column_indices, row in row_indices
		push!(rows, row)
		push!(cols, column)
	end
	return nothing
end

function _bus_input_indices(bus, node_states, incident_states)
	states = node_states[bus]
	return isempty(states) ? incident_states[bus] : states
end

"""
    topology_jacobian_pattern(compiled, components, dev_specs, reference_name)

Build a conservative sparse Jacobian pattern from component-local state blocks and
the electrical network topology. The pattern contains no derivative formulas; it
is intended for colored numerical differentiation of the compiled in-place RHS.
"""
function topology_jacobian_pattern(compiled, components, dev_specs::AbstractDict,
	reference_name::Symbol)
	(; nodes, edges, devs) = components
	nstates = length(unknowns(compiled))

	node_states = Dict{Int, Vector{Int}}()
	for (name, node) in nodes
		bus = extract_v_numbers(name)
		node_states[bus] = _compiled_state_indices(
			compiled, node, jacobian_state_names(Node))
	end

	edge_states = Dict{Symbol, Vector{Int}}()
	device_states = Dict{Symbol, Vector{Int}}()
	incident_states = Dict(bus => Int[] for bus in keys(node_states))

	for (name, edge) in edges
		states = _compiled_state_indices(compiled, edge, jacobian_state_names(Line))
		edge_states[name] = states
		from_bus, to_bus = extract_e_numbers(name)
		append!(incident_states[from_bus], states)
		append!(incident_states[to_bus], states)
	end

	for (name, device) in devs
		spec = dev_specs[name]
		states = _compiled_state_indices(
			compiled, device, jacobian_state_names(spec.model))
		device_states[name] = states
		append!(incident_states[spec.bus], states)
	end

	for states in values(incident_states)
		unique!(states)
	end

	reference_device = devs[reference_name]
	reference_speed = SII.variable_index(compiled, reference_device.ω)
	reference_speed === nothing && error("reference speed is not a compiled state")
	reference_states = [reference_speed]

	rows = Int[]
	cols = Int[]

	for (name, node) in nodes
		bus = extract_v_numbers(name)
		row_indices = node_states[bus]
		column_indices = unique(vcat(
			row_indices,
			incident_states[bus],
			reference_states,
		))
		_append_dense_block!(rows, cols, row_indices, column_indices)
	end

	for (name, edge) in edges
		from_bus, to_bus = extract_e_numbers(name)
		row_indices = edge_states[name]
		column_indices = unique(vcat(
			row_indices,
			_bus_input_indices(from_bus, node_states, incident_states),
			_bus_input_indices(to_bus, node_states, incident_states),
			reference_states,
		))
		_append_dense_block!(rows, cols, row_indices, column_indices)
	end

	for (name, _) in devs
		bus = dev_specs[name].bus
		row_indices = device_states[name]
		column_indices = unique(vcat(
			row_indices,
			_bus_input_indices(bus, node_states, incident_states),
			reference_states,
		))
		_append_dense_block!(rows, cols, row_indices, column_indices)
	end

	# A diagonal entry is retained for every state so implicit solvers can form
	# and factorize W = M - γJ without changing the stored CSC structure.
	append!(rows, 1:nstates)
	append!(cols, 1:nstates)

	pattern = sparse(rows, cols, ones(Float64, length(rows)), nstates, nstates)
	nonzeros(pattern) .= 1.0
	return pattern
end

"""
    check_jacobian_pattern(prob, pattern; atol=1e-8, rtol=1e-7)

Check by central finite differences that the supplied conservative pattern does
not omit a numerically significant Jacobian entry at `prob.u0`.
"""
function check_jacobian_pattern(prob, pattern::SparseMatrixCSC;
	atol::Real = 1e-8, rtol::Real = 1e-7)
	SciMLBase.isinplace(prob) || throw(ArgumentError(
		"check_jacobian_pattern requires an in-place ODEProblem"))
	n = length(prob.u0)
	size(pattern) == (n, n) || throw(DimensionMismatch(
		"Jacobian pattern size $(size(pattern)) does not match $n states"))

	xplus = copy(prob.u0)
	xminus = copy(prob.u0)
	fplus = similar(prob.u0)
	fminus = similar(prob.u0)
	missing = Tuple{Int, Int, Float64}[]
	t = first(prob.tspan)

	for column in 1:n
		x = prob.u0[column]
		h = cbrt(eps(Float64)) * max(1.0, abs(Float64(x)))
		xplus[column] = x + h
		xminus[column] = x - h
		prob.f(fplus, xplus, prob.p, t)
		prob.f(fminus, xminus, prob.p, t)
		column_scale = max(1.0, maximum(abs, fplus), maximum(abs, fminus))
		tolerance = atol + rtol * column_scale
		for row in 1:n
			value = (fplus[row] - fminus[row]) / (2h)
			if abs(value) > tolerance && iszero(pattern[row, column])
				push!(missing, (row, column, Float64(value)))
			end
		end
		xplus[column] = x
		xminus[column] = x
	end

	return (valid = isempty(missing), missing = missing)
end

"""
    prune_jacobian_pattern(prob, pattern; probes=3, perturbation=1e-4,
        atol=1e-11, rtol=1e-9)

Numerically prune a conservative topology pattern.  Only entries already present
in `pattern` are considered, so this cannot introduce nonlocal couplings.  The
union over several deterministic nearby probe points avoids retaining an entire
dense component block merely because a few entries happen to be nonzero at the
operating point.  Diagonal entries are always retained for the implicit solver's
`W = M - γJ` matrix.
"""
function prune_jacobian_pattern(prob, pattern::SparseMatrixCSC;
	probes::Integer = 3, perturbation::Real = 1e-4,
	atol::Real = 1e-11, rtol::Real = 1e-9)
	SciMLBase.isinplace(prob) || throw(ArgumentError(
		"prune_jacobian_pattern requires an in-place ODEProblem"))
	probes >= 1 || throw(ArgumentError("probes must be positive"))
	n = length(prob.u0)
	size(pattern) == (n, n) || throw(DimensionMismatch(
		"Jacobian pattern size $(size(pattern)) does not match $n states"))

	x = copy(prob.u0)
	xplus = copy(x)
	xminus = copy(x)
	fplus = similar(x)
	fminus = similar(x)
	keep = falses(nnz(pattern))
	t = first(prob.tspan)

	for probe in 1:probes
		if probe == 1
			copyto!(x, prob.u0)
		else
			@inbounds for index in eachindex(x)
				scale = max(1.0, abs(Float64(prob.u0[index])))
				x[index] = prob.u0[index] +
					perturbation*scale*sin(index*(probe-1)*0.6180339887498949)
			end
		end

		for column in 1:n
			xvalue = x[column]
			h = cbrt(eps(Float64))*max(1.0, abs(Float64(xvalue)))
			copyto!(xplus, x)
			copyto!(xminus, x)
			xplus[column] = xvalue+h
			xminus[column] = xvalue-h
			prob.f(fplus, xplus, prob.p, t)
			prob.f(fminus, xminus, prob.p, t)
			column_scale = max(1.0, maximum(abs, fplus), maximum(abs, fminus))
			tolerance = atol+rtol*column_scale
			for storage_index in nzrange(pattern, column)
				row = rowvals(pattern)[storage_index]
				value = (fplus[row]-fminus[row])/(2h)
				keep[storage_index] |= abs(value) > tolerance
			end
		end
	end

	rows = Int[]
	cols = Int[]
	for column in 1:n
		for storage_index in nzrange(pattern, column)
			row = rowvals(pattern)[storage_index]
			if keep[storage_index] || row == column
				push!(rows, row)
				push!(cols, column)
			end
		end
	end
	result = sparse(rows, cols, ones(Float64, length(rows)), n, n)
	nonzeros(result) .= 1.0
	return result
end

struct ColoredJacobianCache{F, B, P, Y, J}
	f!::F
	backend::B
	prep::P
	y_cache::Y
	prototype::J
	ncolors::Int
end

"""Prepare a reusable colored numerical Jacobian cache."""
function prepare_colored_jacobian(prob, pattern::SparseMatrixCSC;
	autodiff = AutoFiniteDiff())
	SciMLBase.isinplace(prob) || throw(ArgumentError(
		"prepare_colored_jacobian requires an in-place ODEProblem"))
	size(pattern) == (length(prob.u0), length(prob.u0)) || throw(DimensionMismatch(
		"Jacobian pattern size $(size(pattern)) does not match $(length(prob.u0)) states"))

	structural_pattern = SparseMatrixCSC(
		size(pattern, 1),
		size(pattern, 2),
		copy(pattern.colptr),
		copy(pattern.rowval),
		fill(true, nnz(pattern)),
	)
	backend = AutoSparse(
		autodiff;
		sparsity_detector = ADTypes.KnownJacobianSparsityDetector(structural_pattern),
		coloring_algorithm = SparseMatrixColorings.GreedyColoringAlgorithm(),
	)
	y_cache = similar(prob.u0)
	t0 = first(prob.tspan)
	prep = DI.prepare_jacobian(
		prob.f,
		y_cache,
		backend,
		prob.u0,
		DI.Constant(prob.p),
		DI.Constant(t0),
	)
	colors = SparseMatrixColorings.column_colors(prep)
	prototype = copy(pattern)
	fill!(nonzeros(prototype), 0.0)
	return ColoredJacobianCache(
		prob.f,
		backend,
		prep,
		y_cache,
		prototype,
		isempty(colors) ? 0 : maximum(colors),
	)
end

function colored_jacobian!(J, cache::ColoredJacobianCache, x, p, t)
	DI.jacobian!(
		cache.f!,
		cache.y_cache,
		J,
		cache.prep,
		cache.backend,
		x,
		DI.Constant(p),
		DI.Constant(t),
	)
	return J
end

"""Attach a prepared colored numerical Jacobian to an in-place ODE problem."""
function with_colored_jacobian(prob, cache::ColoredJacobianCache;
	specialization::Type{S} = SciMLBase.NoSpecialize) where {S}
	jac! = let cache = cache
		function (J, x, p, t)
			return colored_jacobian!(J, cache, x, p, t)
		end
	end

	f = SciMLBase.ODEFunction{true, S}(
		prob.f.f;
		mass_matrix = prob.f.mass_matrix,
		jac = jac!,
		jac_prototype = copy(cache.prototype),
		sparsity = cache.prototype,
	)
	return SciMLBase.remake(prob; f)
end
