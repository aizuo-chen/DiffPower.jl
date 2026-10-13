function Line(; name, params::NamedTuple = (;))
	defaults = (
		ωb = 2pi*60,
		rg = 0.015,
		lg = 0.100,
		lg_min = 1e-10,
		k = 1.0,
		id_m = 0.0,
		iq_m = 0.0,
	);
	p = merge(defaults, params)
	lg = max(p.lg, p.lg_min)
	ps = @parameters begin
		ωb = p.ωb
		rg = p.rg
		lg = lg
		k = p.k
	end
	sts = @variables begin
		ωr(t)
		id(t) = p.id_m
		iq(t) = p.iq_m
		isd(t)
		isq(t)
		itd(t)
		itq(t)
		usd(t)
		usq(t)
		utd(t)
		utq(t)
	end
	eqs = [
		D(id) / ωb ~ (usd / k - rg * id - utd) / lg + ωr * iq,
		D(iq) / ωb ~ (usq / k - rg * iq - utq) / lg - ωr * id,
		isd ~ id / k,
		isq ~ iq / k,
		itd ~ id,
		itq ~ iq,
	]
	return System(eqs, t, sts, ps; name)
end

jacobian_state_names(::typeof(Line)) = (:id, :iq)

function Node(; name, params::NamedTuple = (;))
	defaults = (
		ωb = 2π*60,
		g = 0.0,
		g_min = 1e-5,
		b = 0.0,
		ud_m = 0.0,
		uq_m = 1.0,
		cg = 0.0,
	);
	p = merge(defaults, params)
	g = p.g + p.g_min
	ps = @parameters begin
		ωb = p.ωb
		cg = p.cg
		g = g
		b = p.b
	end
	sts = @variables begin
		ωr(t)
		id(t)
		iq(t)
		ud(t) = p.ud_m
		uq(t) = p.uq_m
	end

	if iszero(p.cg)
		y2 = g^2 + (b*ωr)^2
		eqs = [
			ud ~ -(g*id + b*ωr*iq) / y2,
			uq ~ (b*ωr*id - g*iq) / y2,
		]
	else
		eqs = [
			D(ud) / ωb ~ - (id + g*ud - b*uq*ωr) / cg + ωr * uq,
			D(uq) / ωb ~ - (iq + g*uq + b*ud*ωr) / cg - ωr * ud,
		]
	end

	return System(eqs, t, sts, ps; name)
end

jacobian_state_names(::typeof(Node)) = (:ud, :uq)
