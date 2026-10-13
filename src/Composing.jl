function extract_e_numbers(sym::Union{Symbol, String})
	str = string(sym)
	m = match(r"^e(\d+)_(\d+)(?:_(\d+))?$", str)
	if m === nothing
		error("格式错误：必须为 e + 数字_数字 或 e + 数字_数字_序号 形式，如 :e12_13、:e12_13_2")
	end
	return (parse(Int, m[1]), parse(Int, m[2]))
end
function extract_vg_numbers(sym::Union{Symbol, String})
	str = string(sym)
	# 正则解释：
	# ^vg(\d+)   : 以 vg 开头，捕获第一组数字
	# (?:_(\d+))?: 非捕获组，匹配可选的 "_" 后面接数字
	# $          : 字符串结束
	m = match(r"^vg(\d+)(?:_(\d+))?$", str)

	if m === nothing
		error("格式错误：必须为 vg + 数字 或 vg + 数字_数字 形式，例如 :vg12 或 :vg12_1")
	end

	# m[1] 始终是第一组数字
	return parse(Int, m[1])
end
function extract_v_numbers(sym::Union{Symbol, String})
	str = string(sym)
	m = match(r"^v(\d+)$", str)  # 匹配 v + 数字 格式
	if m === nothing
		error("格式错误：必须为 v + 数字 形式，如 :v12")
	end
	return parse(Int, m[1])  # 将字符串转换为整数
end

function composing(vs::Dict{Symbol, System}, es::Dict{Symbol, System}, devs::Dict{Symbol, System}, reference_name::Symbol)
	#网络拓扑连接
	#[:vg1,:vg2,:vg3...,] vg的接口为 ud uq id iq ωr
	#[:v1,:v2,:v3...,] v的接口为 id iq ud uq ωr
	#[:e1_2,:e1_3,:e2_3...] e的接口为 isd isq itd itq usd usq utd utq ωr
	# reference_name 为唯一参考设备名称，例如 :vg1 或 :vg1_1
	# vg_idx 为设备节点编号向量, including GFMs, GFLs, Generators and Infbus.
	# 确保 equations 是一个Equation
	equas = Equation[]
	# 获取参考机的 ODESystem
	haskey(devs, reference_name) || throw(ArgumentError("missing reference device $reference_name"))
	ref_dev = devs[reference_name]
	#  
	n = length(vs)
	id_sums = zeros(Num, n)
	iq_sums = zeros(Num, n)
	# generate all edge-related equations
	for (e, e_sys) in es
		# e 是 Symbol，如 :e1_2, and e_sys 是 e 对象的 ODESystem
		nf, nt = extract_e_numbers(e)
		# 获取 v 的 ODESystem
		vf_sys = vs[Symbol("v", nf)]
		vt_sys = vs[Symbol("v", nt)]
		# 生成连接方程
		push!(equas, e_sys.usd ~ vf_sys.ud)
		push!(equas, e_sys.usq ~ vf_sys.uq)
		push!(equas, e_sys.utd ~ vt_sys.ud)
		push!(equas, e_sys.utq ~ vt_sys.uq)
		push!(equas, e_sys.ωr ~ ref_dev.ω)
		# 计算边对节点电流的净流出
		id_sums[nf] += e_sys.isd
		iq_sums[nf] += e_sys.isq
		id_sums[nt] -= e_sys.itd
		iq_sums[nt] -= e_sys.itq
	end
	# generate all device-realted equations
	for (dev_name, dev) in devs
		ng = extract_vg_numbers(dev_name)
		id_sums[ng] -= dev.id
		iq_sums[ng] -= dev.iq
		# 获取 v 的ODESystem
		v_sys = vs[Symbol("v", ng)]
		push!(equas, dev.ud ~ v_sys.ud)
		push!(equas, dev.uq ~ v_sys.uq)
		if hasproperty(dev, :ωr)
			push!(equas, dev.ωr ~ ref_dev.ω)
		end
	end
	# generate all node-related equations
	for (v, v_sys) in vs
		# v 是 Symbol，例如 :v1,  and v_sys 是 v 对象的 ODESystem
		v_idx = extract_v_numbers(v)
		# 同步频率
		push!(equas, v_sys.ωr ~ ref_dev.ω)
		# 节点电流平衡
		push!(equas, v_sys.id ~ id_sums[v_idx])
		push!(equas, v_sys.iq ~ iq_sums[v_idx])
	end
	return vcat(equas...)
end
