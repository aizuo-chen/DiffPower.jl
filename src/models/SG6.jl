function SG6 end

const SG6_DEFAULTS = SG6_REF_DEFAULTS

function operating_point(::typeof(SG6), voltage, current, params::NamedTuple = (;))
	p = merge(SG6_DEFAULTS, params,
		(id_m=real(current), iq_m=imag(current), ud_m=real(voltage), uq_m=imag(voltage)))
	scale = current_capacity_scale(p.id_m, p.iq_m, p.ratio)
	lg, rg = p.lg/scale, p.rg/scale
	H, Da, Ra = p.H*scale, p.Da*scale, p.Ra/scale
	Xl, Xd, Xq = p.Xl/scale, p.Xd/scale, p.Xq/scale
	Xd′, Xq′, X″ = p.Xd′/scale, p.Xq′/scale, p.X″/scale
	if p.S10 == 0.0 && p.S12 == 0.0
		A, B = 1.0, 0.0
	else
		ρ = sqrt(1.2*p.S12/p.S10)
		A = (ρ-1.2)/(ρ-1.0)
		B = p.S10/(1.0-A)^2
	end
	γqd = (Xq-Xl)/(Xd-Xl)
	vd = real(current)*rg-lg*imag(current)+real(voltage)
	vq = imag(current)*rg+lg*real(current)+imag(voltage)
	ψ″ = hypot(vd-X″*imag(current)+Ra*real(current), vq+X″*real(current)+Ra*imag(current))
	Sat = ψ″ > A ? B*(ψ″-A)^2/ψ″ : 0.0
	Xq_eff = X″+(Xq-X″)/(1+γqd*Sat/ψ″)
	δ = atan(Xq_eff*imag(current)-vd-Ra*real(current), Xq_eff*real(current)+vq+Ra*imag(current))
	rotation = [cos(δ) sin(δ); -sin(δ) cos(δ)]
	Id, Iq = rotation*[real(current); imag(current)]
	Vd, Vq = rotation*[vd; vq]
	ψd, ψq = Vq+Ra*Iq, -Vd-Ra*Id
	ψd″, ψq″ = ψd+X″*Id, ψq+X″*Iq
	ψ″ = hypot(ψd″, ψq″)
	Sat = ψ″ > A ? B*(ψ″-A)^2/ψ″ : 0.0
	Eq′ = ψd″+(Xd′-X″)*Id
	ψd′ = Eq′-(Xd′-Xl)*Id
	Ed′ = -ψq″-(Xq′-X″)*Iq
	ψq′ = Ed′+(Xq′-Xl)*Iq
	Efd = Eq′+(Xd-Xd′)*(Id+(Xd′-X″)/(Xd′-Xl)^2*(Eq′-ψd′-(Xd′-Xl)*Id))+Sat/ψ″*ψd″
	Pmech = ψd*Iq-ψq*Id
	states = (δ=δ, ω=1.0, Eq′=Eq′, ψd′=ψd′, Ed′=Ed′, ψq′=ψq′,
		id=real(current), iq=imag(current))
	parameters = (ωb=p.ωb, Xl=Xl, Xd=Xd, Xq=Xq, Xd′=Xd′, Xq′=Xq′, X″=X″,
		Ra=Ra, Td0′=p.Td0′, Tq0′=p.Tq0′, Td0″=p.Td0″, Tq0″=p.Tq0″,
		H=H, Da=Da, Efd=Efd, Pmech=Pmech, lg=lg, rg=rg, A=A, B=B,
		ratio=p.ratio, scale=scale)
	return (; states, parameters)
end

function SG6(; name, params::NamedTuple = (;))
	p = merge(SG6_DEFAULTS, params)

	(; ωb, lg, rg, H, Da, Ra, Xl, Xd, Xq, Xd′, Xq′, X″, Td0′, Tq0′, Td0″, Tq0″, S10, S12, id_m, iq_m, ud_m, uq_m) = p
	scale = current_capacity_scale(id_m, iq_m, p.ratio)
	lg /= scale
	rg /= scale
	H *= scale
	Da *= scale
	Ra /= scale
	Xl /= scale
	Xd /= scale
	Xq /= scale
	Xd′ /= scale
	Xq′ /= scale
	X″ /= scale

	if S10 == 0.0 && S12 == 0.0
		A, B = 1.0, 0.0
	else
		ρ = sqrt(1.2*S12/S10)
		A = (ρ-1.2)/(ρ-1.0)
		B = S10/(1.0-A)^2
	end

	γqd = (Xq-Xl)/(Xd-Xl)

	vd_m = id_m*rg - lg*iq_m + ud_m
	vq_m = iq_m*rg + lg*id_m + uq_m

	ψ″_m = hypot(vd_m - X″*iq_m + Ra*id_m, vq_m + X″*id_m + Ra*iq_m)
	Sat_m = ψ″_m > A ? B*(ψ″_m-A)^2/ψ″_m : 0.0

	Xq_eff_m = X″ + (Xq-X″)/(1.0 + γqd*Sat_m/ψ″_m)
	δ_m = atan(Xq_eff_m*iq_m - vd_m - Ra*id_m, Xq_eff_m*id_m + vq_m + Ra*iq_m)

	Id_m, Iq_m = [cos(δ_m) sin(δ_m); -sin(δ_m) cos(δ_m)] * [id_m; iq_m]
	Vd_m, Vq_m = [cos(δ_m) sin(δ_m); -sin(δ_m) cos(δ_m)] * [vd_m; vq_m]

	ψd_m = Vq_m + Ra*Iq_m
	ψq_m = -Vd_m - Ra*Id_m

	ψd″_m = ψd_m + X″*Id_m
	ψq″_m = ψq_m + X″*Iq_m
	ψ″_m = hypot(ψd″_m, ψq″_m)
	Sat_m = ψ″_m > A ? B*(ψ″_m-A)^2/ψ″_m : 0.0

	Eq′_m = ψd″_m + (Xd′-X″)*Id_m
	ψd′_m = Eq′_m - (Xd′-Xl)*Id_m

	Ed′_m = -ψq″_m - (Xq′-X″)*Iq_m
	ψq′_m = Ed′_m + (Xq′-Xl)*Iq_m

	LadIfd_m = Eq′_m + (Xd-Xd′)*(Id_m + (Xd′-X″)/(Xd′-Xl)^2*(Eq′_m-ψd′_m-(Xd′-Xl)*Id_m)) + Sat_m/ψ″_m*ψd″_m
	Efd_m = LadIfd_m

	Te_m = ψd_m*Iq_m - ψq_m*Id_m
	Pmech_m = Te_m

	ps = @parameters begin
		ωb = ωb
		Xl = Xl
		Xd = Xd
		Xq = Xq
		Xd′ = Xd′
		Xq′ = Xq′
		X″ = X″
		Ra = Ra
		Td0′ = Td0′
		Tq0′ = Tq0′
		Td0″ = Td0″
		Tq0″ = Tq0″
		H = H
		Da = Da
		Efd = Efd_m
		Pmech = Pmech_m
		lg = lg
		rg = rg
		A = A
		B = B
		ratio = p.ratio
		scale = scale
	end

	sts = @variables begin
		δ(t) = δ_m
		ω(t) = 1.0
		Eq′(t) = Eq′_m
		ψd′(t) = ψd′_m
		Ed′(t) = Ed′_m
		ψq′(t) = ψq′_m

		id(t) = id_m
		iq(t) = iq_m

		ωr(t)
		ud(t)
		uq(t)

		ψd″(t)
		ψq″(t)
		ψd(t)
		ψq(t)
		ψ″(t)
		Sat(t)
		LadIfd(t)

		Te(t)
		pg(t)
		qg(t)

		Id(t)
		Iq(t)

		vd(t)
		vq(t)
		Vd(t)
		Vq(t)
	end

	T1 = [cos(δ) sin(δ); -sin(δ) cos(δ)]
	T2 = [cos(δ) -sin(δ); sin(δ) cos(δ)]

	eq1 = [
		D(δ)/ωb ~ ω - ωr,
		2H*D(ω) ~ (Pmech-Da*(ω-1))/ω - Te,
		Td0′ * D(Eq′) ~ Efd - LadIfd,
		Td0″*D(ψd′) ~ Eq′ - ψd′ - (Xd′-Xl)*Id,
		ψd″ ~ (X″-Xl)/(Xd′-Xl)*Eq′ + (Xd′-X″)/(Xd′-Xl)*ψd′,
		Tq0′ * D(Ed′) ~ -Ed′ + (Xq-Xq′)*(Iq - (Xq′-X″)/(Xq′-Xl)^2*(Ed′-ψq′+(Xq′-Xl)*Iq)) + (Xq-Xl)/(Xd-Xl)*Sat/ψ″*ψq″,
		Tq0″*D(ψq′) ~ Ed′ - ψq′ + (Xq′-Xl)*Iq,
		ψq″ ~ -((X″-Xl)/(Xq′-Xl)*Ed′ + (Xq′-X″)/(Xq′-Xl)*ψq′),
		ψ″ ~ sqrt(ψd″^2 + ψq″^2),
		Sat ~ ifelse(ψ″ > A, B*(ψ″-A)^2/ψ″, 0.0),
		LadIfd ~ Eq′ + (Xd-Xd′)*(Id + (Xd′-X″)/(Xd′-Xl)^2*(Eq′-ψd′-(Xd′-Xl)*Id)) + Sat/ψ″*ψd″, ψd ~ ψd″ - X″*Id,
		ψq ~ ψq″ - X″*Iq,
		Te ~ ψd*Iq - ψq*Id,
		Vd ~ -ω*ψq - Ra*Id,
		Vq ~ ω*ψd - Ra*Iq,
		pg ~ Vd*Id + Vq*Iq,
		qg ~ -Vd*Iq + Vq*Id,
		D(id) / ωb ~ -rg*id/lg + (vd-ud)/lg + ωr*iq,
		D(iq)/ωb ~ -rg*iq/lg + (vq-uq)/lg - ωr*id,
	]

	eq2 = [Id; Iq] .~ T1*[id; iq]
	eq3 = [vd; vq] .~ T2*[Vd; Vq]

	eqs = vcat(eq1, eq2, eq3)

	return System(eqs, t, sts, ps; name = name)
end
