function calc_Ybus(
	bus::Matrix{Float64}, branch::Matrix{Float64};
	baseMVA::Float64 = 100.0,
	g_min::Float64 = 1e-5,
	lg_min::Float64 = 1e-10,
)::Matrix{ComplexF64}
	n = size(bus, 1)                                 # 母线数量
	Y = zeros(ComplexF64, n, n)                      # 稠密导纳矩阵

	# 1. 处理支路（branch）导纳
	@inbounds for k in axes(branch, 1)
		i   = Int(branch[k, 1])        # from
		j   = Int(branch[k, 2])        # to
		r   = branch[k, 3]
		x   = max(branch[k, 4], lg_min)
		bc  = branch[k, 5]             # 总充电电纳（整条线路）
		tap = branch[k, 9]             # 变比，0 表示 1.0
		# shift = branch[k, 10] * π/180    # 移相角（度→弧度）
		shift = 0.0

		tap = tap == 0 ? 1.0 : tap
		a   = tap * exp(1im * shift)      # 复数变比

		ys = 1 / (r + 1im*x)              # 串联导纳
		yc = 1im * bc / 2                 # 每端充电导纳（π型）

		Yii = (ys + yc) / (tap^2)
		Yjj = ys + yc
		Yij = -ys / conj(a)
		Yji = -ys / a

		Y[i, i] += Yii
		Y[j, j] += Yjj
		Y[i, j] += Yij
		Y[j, i] += Yji
	end

	# 2. 添加母线并联导纳 (shunt)
	@inbounds for i in 1:n
		Gs = bus[i, 5]   # 单位：MW
		Bs = bus[i, 6]   # 单位：MVAR
		Y[i, i] += (Gs + 1im * Bs) / baseMVA + g_min   # 转换为 p.u. 导纳
	end

	return Y
end

function calc_cap(branch::Matrix{Float64})::Vector{Float64}
	n = Int(maximum(branch[:, 1:2]))
	Bc = zeros(Float64, n)  # 列向量，长度 n

	# 遍历每条支路，将 bc/2 加到 from 和 to 节点
	@inbounds for k in axes(branch, 1)
		i = Int(round(branch[k, 1]))
		j = Int(round(branch[k, 2]))
		bc = branch[k, 5]       # 总充电电纳（p.u.）
		tap = branch[k, 9] == 0 ? 1.0 : branch[k, 9]
		# Match calc_Ybus: the from-side shunt is referred through tap².
		Bc[i] += bc / (2.0*tap^2)
		Bc[j] += bc / 2.0
	end

	return Bc  # 返回列向量，单位：p.u.（容性电纳）
end

function calc_injections(bus::Matrix{Float64}, gen::Matrix{Float64}; baseMVA::Float64 = 100.0)::Tuple{Vector{Float64}, Vector{Float64}}
	P_inj = -bus[:, 3] ./ baseMVA   # 负荷有功
	Q_inj = -bus[:, 4] ./ baseMVA   # 负荷无功

	@inbounds for g in axes(gen, 1)
		i = Int(gen[g, 1])          # 发电机所在母线编号
		P_inj[i] += gen[g, 2] / baseMVA
	end
	return P_inj, Q_inj
end

struct PFParams
	Ybus::Matrix{ComplexF64}
	P_inj::Vector{Float64}
	Q_inj::Vector{Float64}
	bus_type::Vector{Int}
	V_set::Vector{Float64}
end

function residual(x::AbstractVector{T}, p::PFParams) where T
	n_bus = length(p.P_inj)
	θ = @view x[1:n_bus]
	V = @view x[(n_bus+1):(2n_bus)]
	U = V .* cis.(θ)
	S = U .* conj.(p.Ybus * U)
	res = zeros(T, 2n_bus)

	@inbounds for i in 1:n_bus
		if p.bus_type[i] == 3       # slack: angle and voltage magnitude
			res[i] = θ[i]
			res[n_bus+i] = V[i] - p.V_set[i]
		elseif p.bus_type[i] == 2   # PV: active power and voltage magnitude
			res[i] = real(S[i]) - p.P_inj[i]
			res[n_bus+i] = V[i] - p.V_set[i]
		else                        # PQ: active and reactive power
			res[i] = real(S[i]) - p.P_inj[i]
			res[n_bus+i] = imag(S[i]) - p.Q_inj[i]
		end
	end

	return res
end

function newton_pf(
	bus::Matrix{Float64}, branch::Matrix{Float64}, gen::Matrix{Float64};
	tol::Float64 = 1e-10,
	maxiter::Int = 20,
	baseMVA::Float64 = 100.0,
	g_min::Float64 = 1e-5,
	lg_min::Float64 = 1e-10,
)

	slack_bus = findfirst(bus[:, 2] .== 3) # Slack节点
	isnothing(slack_bus) && error("power flow requires one slack bus (bus type 3)")
	# 构建节点导纳矩阵
	Ybus = calc_Ybus(bus, branch; baseMVA = baseMVA, g_min = g_min, lg_min = lg_min)
	# 计算注入功率
	P_inj, Q_inj = calc_injections(bus, gen; baseMVA = baseMVA)
	params = PFParams(Ybus, P_inj, Q_inj, Int.(bus[:, 2]), copy(bus[:, 8]))
	θ = deg2rad.(bus[:, 9] .- bus[slack_bus, 9])
	V = copy(bus[:, 8])
	x = [θ; V]

	for iter ∈ 1:maxiter
		f = residual(x, params)
		fnorm = maximum(abs.(f))

		println("iter $iter, max residual = $fnorm")
		if fnorm < tol
			println("Converged!")
			θ = @view x[1:size(bus, 1)]
			V = @view x[(size(bus, 1)+1):(2size(bus, 1))]
			# Return standard power-flow phasors. The EMT reference-axis rotation is
			# applied once, later, by align_reference! for the selected reference model.
			V_bus = V .* cis.(θ)
			# === 计算支路电流 ===
			from_idx = Int.(branch[:, 1])
			to_idx = Int.(branch[:, 2])
			V_f = V_bus[from_idx]
			V_t = V_bus[to_idx]
			x_branch = max.(branch[:, 4], lg_min)
			ys = 1.0 ./ (branch[:, 3] .+ 1im .* x_branch)
			tap = ifelse.(branch[:, 9] .== 0, 1.0, branch[:, 9])
			I_branch = ys .* (V_f ./ tap .- V_t)
			I_bus = Ybus*V_bus
			Pds = bus[:, 3]/baseMVA
			Qds = bus[:, 4]/baseMVA
			# The MATPOWER bus shunt (Gs + jBs) is a fixed nodal admittance.
			# Only the charging susceptance at the two ends of each line is
			# represented by Node's dynamic capacitance.
			Yds = (Pds-1im*Qds) ./ abs.(V_bus) .^ 2 .+
				(bus[:, 5] .+ 1im .* bus[:, 6]) ./ baseMVA
			I_gen = I_bus + conj((bus[:, 3]+1im*bus[:, 4])/baseMVA ./ V_bus)
			c_bus = calc_cap(branch)
			return (
				V_bus = V_bus,
				I_gen = I_gen,
				Yds = Yds,
				I_branch = I_branch,
				c_bus = c_bus,
				slack_bus = slack_bus,
			)
		end

		J = ForwardDiff.jacobian(x -> residual(x, params), x)
		Δx = -J \ f
		x += Δx
	end
	error("power flow did not converge in $maxiter iterations")
end
