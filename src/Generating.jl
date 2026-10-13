device_params(spec) = get(spec, :params, (;))

function validate_device_spec(name, spec, bus_count)
	hasproperty(spec, :bus) || throw(ArgumentError("device $name is missing `bus`"))
	hasproperty(spec, :model) || throw(ArgumentError("device $name is missing `model`"))
	hasproperty(spec, :power) || throw(ArgumentError("device $name is missing `power`"))
	bus = spec.bus
	1 <= bus <= bus_count || throw(ArgumentError("device $name uses invalid bus $bus"))
	extract_vg_numbers(name) == bus || throw(ArgumentError(
		"device name $name maps to bus $(extract_vg_numbers(name)), but its spec uses bus $bus"))
	return nothing
end

function resolve_device_powers(pf, dev_specs::AbstractDict)
	groups = Dict{Int, Vector{Symbol}}()
	for (name, spec) in dev_specs
		validate_device_spec(name, spec, length(pf.V_bus))
		push!(get!(groups, spec.bus, Symbol[]), name)
	end
	powers = Dict{Symbol, ComplexF64}()
	for (bus, names) in groups
		remainders = filter(name -> isnothing(dev_specs[name].power), names)
		length(remainders) == 1 || throw(ArgumentError(
			"bus $bus must have exactly one device with `power=nothing`; found $(length(remainders))"))
		specified = zero(ComplexF64)
		for name in names
			power = dev_specs[name].power
			isnothing(power) && continue
			power isa Number || throw(ArgumentError("device $name power must be a number or nothing"))
			powers[name] = ComplexF64(power)
			specified += power
		end
		total = pf.V_bus[bus] * conj(pf.I_gen[bus])
		powers[only(remainders)] = total - specified
	end
	return powers
end

device_current(pf, name, dev_specs, powers) =
	conj(powers[name] / pf.V_bus[dev_specs[name].bus])

function align_reference!(pf, dev_specs::AbstractDict, reference_name::Symbol)
	haskey(dev_specs, reference_name) || throw(ArgumentError("missing reference device $reference_name"))
	reference = dev_specs[reference_name]
	reference.bus == pf.slack_bus || throw(ArgumentError(
		"reference device $reference_name is on bus $(reference.bus), not slack bus $(pf.slack_bus)"))
	powers = resolve_device_powers(pf, dev_specs)
	current = device_current(pf, reference_name, dev_specs, powers)
	reference_phasor = reference_axis_phasor(
		reference.model, pf.V_bus[reference.bus], current, device_params(reference))
	rotation = cis(π/2 - angle(reference_phasor))
	pf.V_bus .*= rotation
	pf.I_gen .*= rotation
	pf.I_branch .*= rotation
	return nothing
end

function generate_edges(branch, I_branch; lg_min = 1e-10)
	edges = Dict{Symbol, System}()
	params_dict = Dict{Symbol, NamedTuple}()
	endpoint_counts = Dict{Tuple{Int, Int}, Int}()
	for idx in axes(branch, 1)
		node1, node2 = Int(branch[idx, 1]), Int(branch[idx, 2])
		endpoints = (node1, node2)
		parallel_index = get(endpoint_counts, endpoints, 0) + 1
		endpoint_counts[endpoints] = parallel_index
		name = parallel_index == 1 ?
			Symbol("e$(node1)_$(node2)") :
			Symbol("e$(node1)_$(node2)_$(parallel_index)")
		params = (id_m = real(I_branch[idx]), iq_m = imag(I_branch[idx]),
			rg = branch[idx, 3], lg = branch[idx, 4], lg_min,
			k = branch[idx, 9] == 0 ? 1 : branch[idx, 9])
		params_dict[name] = params
		edges[name] = Line(; name, params)
	end
	return params_dict, edges
end

function generate_nodes(c_bus, Yds, V_bus; g_min = 1e-5)
	nodes = Dict{Symbol, System}()
	params_dict = Dict{Symbol, NamedTuple}()
	for idx in eachindex(V_bus)
		isfinite(c_bus[idx]) && c_bus[idx] >= 0 || throw(ArgumentError(
			"bus $idx has invalid line-charging capacitance $(c_bus[idx])"))
		name = Symbol("v$(idx)")
		params = (g = real(Yds[idx]), g_min, b = imag(Yds[idx]),
			ud_m = real(V_bus[idx]), uq_m = imag(V_bus[idx]), cg = c_bus[idx])
		params_dict[name] = params
		nodes[name] = Node(; name, params)
	end
	return params_dict, nodes
end

function generate_devs(pf, dev_specs::AbstractDict; params::AbstractDict = Dict())
	devs = Dict{Symbol, System}()
	params_dict = Dict{Symbol, NamedTuple}()
	powers = resolve_device_powers(pf, dev_specs)
	currents = Dict{Symbol, ComplexF64}()
	for (name, spec) in dev_specs
		current = device_current(pf, name, dev_specs, powers)
		currents[name] = current
		model_params = merge(device_params(spec), get(params, name, (;)))
		initial_params = merge((id_m = real(current), iq_m = imag(current),
				ud_m = real(pf.V_bus[spec.bus]), uq_m = imag(pf.V_bus[spec.bus])), model_params)
		params_dict[name] = initial_params
		devs[name] = spec.model(; name, params = initial_params)
	end
	return params_dict, devs, powers, currents
end

function generate_components(pf, branch, dev_specs::AbstractDict;
	params::AbstractDict = Dict(), g_min = 1e-5, lg_min = 1e-10)
	dev_params, devs, dev_powers, dev_currents = generate_devs(pf, dev_specs; params)
	node_params, nodes = generate_nodes(pf.c_bus, pf.Yds, pf.V_bus; g_min)
	edge_params, edges = generate_edges(branch, pf.I_branch; lg_min)
	return (; devs, nodes, edges, dev_params, dev_powers, dev_currents,
		node_params, edge_params)
end

function remake_operating_point(prob, devs::AbstractDict, pf, dev_specs::AbstractDict;
	params::AbstractDict = Dict(), names = keys(dev_specs))
	powers = resolve_device_powers(pf, dev_specs)
	state_pairs = Pair[]
	parameter_pairs = Pair[]
	for name in names
		haskey(devs, name) || throw(ArgumentError("missing generated model $name"))
		haskey(dev_specs, name) || throw(ArgumentError("missing device specification $name"))
		spec = dev_specs[name]
		current = device_current(pf, name, dev_specs, powers)
		model_params = merge(device_params(spec), get(params, name, (;)))
		point = operating_point(spec.model, pf.V_bus[spec.bus], current, model_params)
		mapped = operating_point_pairs(devs[name], point)
		append!(state_pairs, mapped.states)
		append!(parameter_pairs, mapped.parameters)
	end
	return SciMLBase.remake(prob; u0 = state_pairs, p = parameter_pairs)
end
