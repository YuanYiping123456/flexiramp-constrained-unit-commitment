function calculate_flexiramp_capacity(units, x₀, p₀, NS, NG, NT)
	flexiramp_up, flexiramp_dw = zeros(NS, NT), zeros(NS, NT)
	for s in 1:NS
		tem1, tem2 = zeros(NG, NT), zeros(NS, NT)
		for i in 1:NG
			for t in 1:NT
				tem1[i, t] = units.p_max[i, 1] * x₀[i, t] - p₀[(s - 1) * NG + i, t]
				tem2[i, t] = p₀[(s - 1) * NG + i, t] - units.p_min[i, 1] * x₀[i, t]
			end
		end
		flexiramp_up[s, :] = sum(tem1, dims = 1)
		flexiramp_dw[s, :] = sum(tem2, dims = 1)
	end
	return flexiramp_up, flexiramp_dw
end
