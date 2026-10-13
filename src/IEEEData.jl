function make_continuous_bus(bus::Matrix, branch::Matrix, gen::Matrix)
	# Collect all bus ids from bus matrix
	all_bus_ids = Set{Int}()
	for i in 1:size(bus, 1)
		push!(all_bus_ids, Int(bus[i, 1]))
	end

	# Collect buses that appear in branches
	connected_bus_ids = Set{Int}()
	if size(branch, 1) > 0
		for i in 1:size(branch, 1)
			push!(connected_bus_ids, Int(branch[i, 1]))  # from bus
			push!(connected_bus_ids, Int(branch[i, 2]))  # to bus
		end
	end

	# Check for isolated buses
	isolated = setdiff(all_bus_ids, connected_bus_ids)
	if !isempty(isolated)
		error("Isolated buses detected! Bus numbers: $(sort(collect(isolated))).\n" *
			  "No isolated buses allowed. Please fix the network topology first.")
	end

	# Create mapping: old bus id → new continuous id (1 to nbus)
	sorted_bus_ids = sort(collect(all_bus_ids))
	nbus = length(sorted_bus_ids)
	old_to_new = Dict{Int, Int}(old => new for (new, old) in enumerate(sorted_bus_ids))

	# Create new bus matrix with continuous numbering
	new_bus = copy(bus)
	for i in axes(new_bus, 1)
		new_bus[i, 1] = old_to_new[Int(new_bus[i, 1])]
	end

	# Update branch matrix (f_bus and t_bus)
	new_branch = copy(branch)
	if size(new_branch, 1) > 0
		for i in axes(new_branch, 1)
			new_branch[i, 1] = old_to_new[Int(new_branch[i, 1])]
			new_branch[i, 2] = old_to_new[Int(new_branch[i, 2])]
		end
	end

	# Update gen matrix (GEN_BUS column)
	new_gen = copy(gen)
	if size(new_gen, 1) > 0
		for i in axes(new_gen, 1)
			new_gen[i, 1] = old_to_new[Int(new_gen[i, 1])]
		end
	end

	println("Success: Bus renumbering completed! No isolated buses. " *
			"Buses now numbered continuously from 1 to $nbus.")

	return new_bus, new_branch, new_gen, old_to_new
end

function make_continuous_bus!(bus::Matrix, branch::Matrix, gen::Matrix)
	# Collect all bus ids from bus matrix
	all_bus_ids = Set{Int}()
	for i in 1:size(bus, 1)
		push!(all_bus_ids, Int(bus[i, 1]))
	end

	# Collect buses that appear in branches
	connected_bus_ids = Set{Int}()
	if size(branch, 1) > 0
		for i in 1:size(branch, 1)
			push!(connected_bus_ids, Int(branch[i, 1]))  # from bus
			push!(connected_bus_ids, Int(branch[i, 2]))  # to bus
		end
	end

	# Check for isolated buses
	isolated = setdiff(all_bus_ids, connected_bus_ids)
	if !isempty(isolated)
		error("Isolated buses detected! Bus numbers: $(sort(collect(isolated))).")
	end

	# Create mapping: old bus id → new continuous id (1 to nbus)
	sorted_bus_ids = sort(collect(all_bus_ids))
	old_to_new = Dict{Int, Int}(old => new for (new, old) in enumerate(sorted_bus_ids))
	nbus = length(sorted_bus_ids)

	# ======================
	# 直接原地修改！无 copy！
	# ======================
	for i in axes(bus, 1)
		bus[i, 1] = old_to_new[Int(bus[i, 1])]
	end

	if size(branch, 1) > 0
		for i in axes(branch, 1)
			branch[i, 1] = old_to_new[Int(branch[i, 1])]
			branch[i, 2] = old_to_new[Int(branch[i, 2])]
		end
	end

	if size(gen, 1) > 0
		for i in axes(gen, 1)
			gen[i, 1] = old_to_new[Int(gen[i, 1])]
		end
	end

	println("Success: In-place bus renumbering done! Buses 1→$nbus.")
	return old_to_new  # 只返回映射，不返回矩阵
end
