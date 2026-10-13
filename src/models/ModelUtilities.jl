## General design of filter parameters:
# lf = 1/cf/fc^2;
# rf = lf/3;
function current_capacity_scale(id, iq, ratio)
	ratio > 0 || throw(ArgumentError("ratio must be positive"))
	current_magnitude = hypot(id, iq)
	return current_magnitude > 1 ? ratio * ceil(current_magnitude) : 1.0
end

function operating_point end

function operating_point_pairs(system, point)
	state_pairs = Pair[]
	parameter_pairs = Pair[]
	for (name, value) in pairs(point.states)
		push!(state_pairs, getproperty(system, name) => value)
	end
	for (name, value) in pairs(point.parameters)
		push!(parameter_pairs, getproperty(system, name) => value)
	end
	return (; states = state_pairs, parameters = parameter_pairs)
end
