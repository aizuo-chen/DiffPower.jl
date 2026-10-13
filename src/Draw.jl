function dev_trajs(eles, sols; devs, edges)
	results = Dict{Symbol, Any}()
	results[:t] = sols.t
	for ele in eles
		if string(ele)[1]=='v'
			U1 = sqrt.(sols[devs[ele].ucd] .^ 2 + sols[devs[ele].ucq] .^ 2);
			I1 = sqrt.(sols[devs[ele].i1d] .^ 2 + sols[devs[ele].i1q] .^ 2);
			Pg = sols[devs[ele].pg];
			Qg = sols[devs[ele].qg];
			ω = sols[devs[ele].ω];
			results[ele] = (U1 = U1, I1 = I1, Pg = Pg, Qg = Qg, ω = ω)
		else
			Is = sqrt.(sols[edges[ele].isd] .^ 2 + sols[edges[ele].isq] .^ 2);
			It = sqrt.(sols[edges[ele].itd] .^ 2 + sols[edges[ele].itq] .^ 2);
			results[ele] = (Is = Is, It = It)
		end
	end
	results
end
