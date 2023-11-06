function calculate_netload_curve(ben_pᵩ, pᵩ, re_pᵩ, ben_pᵨ, pᵨ, re_pᵨ, loads, winds)
	netload_curve_1, netload_curve_2, netload_curve_3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	tem1, tem2, tem3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	str1, str2, str3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	for s in 1:NS
		t1 = (s - 1) * NW + 1
		t2 = s * NW
		tem1[s, :] .= sum(ben_pᵩ[t1:t2, :], dims = 1)[1, :]
		tem2[s, :] .= sum(pᵩ[t1:t2, :], dims = 1)[1, :]
		tem3[s, :] .= sum(re_pᵩ[t1:t2, :], dims = 1)[1, :]
		str1[s, :] .= sum(pᵨ[t1:t2, :], dims = 1)[1, :]
		str2[s, :] .= sum(re_pᵨ[t1:t2, :], dims = 1)[1, :]
		str3[s, :] .= sum(ben_pᵨ[t1:t2, :], dims = 1)[1, :]
	end

	for s in 1:NS
		netload_curve_1[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem1[s, :] - str1[s, :]
		netload_curve_2[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem2[s, :] - str2[s, :]
		netload_curve_3[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[s, :] .* sum(winds.p_max) + tem3[s, :] - str3[s, :]
	end

	return netload_curve_1, netload_curve_2, netload_curve_3
end

function calculate_ini_netload_curve(ben_pᵩ, pᵩ, re_pᵩ, ben_pᵨ, pᵨ, re_pᵨ, loads, winds)
	netload_curve_1, netload_curve_2, netload_curve_3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	tem1, tem2, tem3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	str1, str2, str3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
	for s in 1:NS
		t1 = (s - 1) * NW + 1
		t2 = s * NW
		tem1[s, :] .= sum(ben_pᵩ[t1:t2, :], dims = 1)[1, :]
		tem2[s, :] .= sum(pᵩ[t1:t2, :], dims = 1)[1, :]
		tem3[s, :] .= sum(re_pᵩ[t1:t2, :], dims = 1)[1, :]
		str1[s, :] .= sum(pᵨ[t1:t2, :], dims = 1)[1, :]
		str2[s, :] .= sum(re_pᵨ[t1:t2, :], dims = 1)[1, :]
		str3[s, :] .= sum(ben_pᵨ[t1:t2, :], dims = 1)[1, :]
	end

	for s in 1:NS
		netload_curve_1[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max)
		netload_curve_2[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max)
		netload_curve_3[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[s, :] .* sum(winds.p_max)
	end

	return netload_curve_1, netload_curve_2, netload_curve_3
end

# flexiramp_up_requirement, flexiramp_dw_requirement = calculate_UpandDw_flexramp_requirement(netload_curve, NT, NS)

function calculate_UpandDw_flexramp_requirement(netload_curve, NT, NS)
	flexiramp_up_requirement, flexiramp_dw_requirement = zeros(NS, NT), zeros(NS, NT)
	for t in 2:NT
		for s in 1:NS
			println([t s])
			flexiramp_up_requirement[s, t - 1] = max(0, netload_curve[s, t] - netload_curve[s, t - 1])
			flexiramp_dw_requirement[s, t - 1] = abs(min(0, netload_curve[s, t] - netload_curve[s, t - 1]))
		end
	end
	return flexiramp_up_requirement, flexiramp_dw_requirement
end
