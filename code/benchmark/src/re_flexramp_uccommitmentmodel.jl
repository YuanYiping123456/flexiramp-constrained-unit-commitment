using JuMP, Gurobi, Test, DelimitedFiles

include("linearization.jl")
include("powerflowcalculation.jl")

function re_flexramp_scucmodel(NT::Int64, NB::Int64, NG::Int64, ND::Int64, NC::Int64, units::unit, loads::load, winds::wind, lines::transmissionline, config_param::config)
	println("Step-3: Creating dispatching model")

	if config_param.is_NetWorkCon == 1
		Adjacmatrix_BtoG, Adjacmatrix_B2D, Gsdf = linearpowerflow(units, lines, loads, NG, NB, ND, NL)
		Adjacmatrix_BtoW = zeros(NB, length(winds.index))
		for i in 1:length(winds.index)
			Adjacmatrix_BtoW[winds.index[i, 1], i] = 1
		end
	end

	NS = winds.scenarios_nums
	NW = length(winds.index)

	# creat scucsimulation_model
	scuc = Model(Gurobi.Optimizer)

	# NS = 1 # for testModel(Gurobi.Optimizer)

	# binary variables
	@variable(scuc, x[1:NG, 1:NT], Bin)
	@variable(scuc, u[1:NG, 1:NT], Bin)
	@variable(scuc, v[1:NG, 1:NT], Bin)

	# continuous variables
	@variable(scuc, pg₀[1:(NG * NS), 1:NT]>=0)
	@variable(scuc, pgₖ[1:(NG * NS), 1:NT, 1:3]>=0)
	@variable(scuc, su₀[1:NG, 1:NT]>=0)
	@variable(scuc, sd₀[1:NG, 1:NT]>=0)
	@variable(scuc, sr⁺[1:(NG * NS), 1:NT]>=0)
	@variable(scuc, sr⁻[1:(NG * NS), 1:NT]>=0)
	@variable(scuc, Δpd[1:(ND * NS), 1:NT]>=0)
	@variable(scuc, Δpw[1:(NW * NS), 1:NT]>=0)

	# # pss variables
	# @variable(scuc, κ⁺[1:(NC * NS), 1:NT], Bin) # charge status
	# @variable(scuc, κ⁻[1:(NC * NS), 1:NT], Bin) # discharge status
	# @variable(scuc, pc⁺[1:(NC * NS), 1:NT] >= 0)# charge power
	# @variable(scuc, pc⁻[1:(NC * NS), 1:NT] >= 0)# discharge power
	# @variable(scuc, qc[1:(NC * NS), 1:NT] >= 0) # cumsum power
	# # @variable(scuc, pss_sumchargeenergy[1:NC * NS, 1] >= 0)

	# @variable(scuc, ζ[1:NG, 1:NT] >= 0)
	# @variable(scuc, z[1:NG^2, 1:NT], Bin)

	refcost, eachslope = linearizationfuelcurve(units, NG)

	c₀ = config_param.is_CoalPrice
	pₛ = scenarios_prob
	plentycoffi_1 = config_param.is_LoadsCuttingCoefficient * 1e10
	plentycoffi_2 = config_param.is_WindsCuttingCoefficient
	ρ⁺ = c₀ * 2
	ρ⁻ = c₀ * 2

	println("start...")
	println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

	# model-1: MIQP with quartic fuel equation
	# @objective(scuc, Min, sum(sum(su₀[i, t] + sd₀[i, t] for i in 1:NG) for t in 1:NT)+
	#     pₛ*c₀/100*(sum(sum(sum(sum(units.coffi_a[i,1] * (pg₀[i + (s - 1) * NG,t]^2) + units.coffi_b[i,1] * pg₀[i + (s - 1) * NG,t] + units.coffi_c[i,1] for t in 1:NT)) for s in 1:NS) for i in 1:NG))+
	#     pₛ*plentycoffi_1*sum(sum(sum(Δpd[1+(s-1)*ND : s*ND, t]) for t in 1:NT) for s in 1:NS)+
	#     pₛ*plentycoffi_2*sum(sum(sum(Δpw[1+(s-1)*NW : s*NW, t]) for t in 1:NT) for s in 1:NS))

	# model-2:MILP with piece linearization equation of nonliear equation
	@objective(scuc,
		Min,
		100*sum(sum(su₀[i, t] + sd₀[i, t] for i in 1:NG) for t in 1:NT)+
		pₛ*
		c₀*
		(sum(sum(sum(sum(pgₖ[i + (s - 1) * NG, t, :] .* eachslope[:, i] for t in 1:NT))
				 for s in 1:NS) for i in 1:NG)+
		sum(sum(sum(x[:, t] .* refcost[:, 1] for t in 1:NT)) for s in 1:NS))+
		pₛ*
		plentycoffi_1*
		sum(sum(sum(Δpd[(1 + (s - 1) * ND):(s * ND), t]) for t in 1:NT) for s in 1:NS)+
		pₛ*
		plentycoffi_2*
		sum(sum(sum(Δpw[(1 + (s - 1) * NW):(s * NW), t]) for t in 1:NT) for s in 1:NS))

	#
	# for test
	# @objective(scuc, Min, 0)
	println("objective_function")
	println("\t MILP_type objective_function \t\t\t\t\t\t done")

	println("subject to.")

	# minimum shutup and shutdown ductration limits
	onoffinit, Lupmin, Ldownmin = zeros(NG, 1), zeros(NG, 1), zeros(NG, 1)
	for i in 1:NG
		onoffinit[i] = ((units.x_0[i, 1] > 0.5) ? 1 : 0)
		Lupmin[i] = min(NT, units.min_shutup_time[i] * onoffinit[i])
		Ldownmin[i] = min(NT, (units.min_shutdown_time[i, 1]) * (1 - onoffinit[i]))
	end

	@constraint(scuc,
		[i = 1:NG, t = 1:Int64((Lupmin[i] + Ldownmin[i]))],
		x[i, t]==onoffinit[i])

	for i in 1:NG
		for t in Int64(max(1, Lupmin[i])):NT
			LB = Int64(max(t - units.min_shutup_time[i, 1] + 1, 1))
			@constraint(scuc, sum(u[i, r] for r in LB:t)<=x[i, t])
		end
		for t in Int64(max(1, Ldownmin[i])):NT
			LB = Int64(max(t - units.min_shutup_time[i, 1] + 1, 1))
			@constraint(scuc, sum(v[i, r] for r in LB:t)<=(1 - x[i, t]))
		end
	end
	println("\t constraints: 1) minimum shutup/shutdown time limits\t\t\t done")

	# binary variable logic
	@constraint(scuc,
		[i = 1:NG, t = 1:NT],
		u[i, t] - v[i, t]==x[i, t] - ((t == 1) ? onoffinit[i] : x[i, t - 1]))
	@constraint(scuc, [i = 1:NG, t = 1:NT], u[i, t] + v[i, t]<=1)
	println("\t constraints: 2) binary variable logic\t\t\t\t\t done")

	# shutup/shutdown cost
	shutupcost = units.coffi_cold_shutup_1
	shutdowncost = units.coffi_cold_shutdown_1
	@constraint(scuc, [t = 1], su₀[:, t].>=shutupcost .* (x[:, t] - onoffinit[:, 1]))
	@constraint(scuc, [t = 1], sd₀[:, t].>=shutdowncost .* (onoffinit[:, 1] - x[:, t]))
	@constraint(scuc, [t = 2:NT], su₀[:, t].>=shutupcost .* u[:, t])
	@constraint(scuc, [t = 2:NT], sd₀[:, t].>=shutdowncost .* v[:, t])
	println("\t constraints: 3) shutup/shutdown cost\t\t\t\t\t done")

	# loadcurtailments and spoliedwinds limits
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		Δpw[(1 + (s - 1) * NW):(s * NW), t].<=
		winds.scenarios_curve[s, t] * winds.p_max[:, 1])
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		Δpd[(1 + (s - 1) * ND):(s * ND), t].<=loads.load_curve[:, t])
	# @constraint(scuc, [s=1:NS, t = 1:NT], Δpw[1+(s-1)*NW:s*NW, t] .== zeros(NW,1))
	# @constraint(scuc, [s=1:NS, t = 1:NT], Δpd[1+(s-1)*ND:s*ND, t] .== zeros(ND,1))
	println("\t constraints: 4) loadcurtailments and spoliedwinds\t\t\t done")

	# generatos power limits
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		pg₀[(1 + (s - 1) * NG):(s * NG), t] + sr⁺[(1 + (s - 1) * NG):(s * NG), t].==
		units.p_max[:, 1] .* x[:, t])
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		pg₀[(1 + (s - 1) * NG):(s * NG), t] - sr⁻[(1 + (s - 1) * NG):(s * NG), t].==
		units.p_min[:, 1] .* x[:, t])
	println("\t constraints: 5) generatos power limits\t\t\t\t\t done")

	# system reserves
	# forcast_error = 0.05
	# forcast_reserve = winds.scenarios_curve * sum(winds.p_max[:,1]) * forcast_error
	# @constraint(scuc,[s = 1:NS, t = 1:NT, i = 1:NG], sum(sr⁺[1 + (s - 1) * NG:s * NG, t]) >= 0.05 * units.p_max[i,1] * x[i,t])
	# @constraint(scuc,[s = 1:NS, t = 1:NT], sum(sr⁻[1 + (s - 1) * NG:s * NG, t]) >= config_param.is_Alpha * forcast_reserve[s,t] + config_param.is_Belta * sum(loads.load_curve[:,t]))

	forcast_error = 0.05
	forcast_reserve = winds.scenarios_curve * sum(winds.p_max[:, 1]) * forcast_error
	@constraint(scuc,
		[s = 1:NS, t = 1:NT, i = 1:NG],
		sum(sr⁺[(1 + (s - 1) * NG):(s * NG), t])>=0.005 * units.p_max[i, 1] * x[i, t])
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		sum(sr⁻[(1 + (s - 1) * NG):(s * NG), t])>=
		0.005 * (config_param.is_Alpha * forcast_reserve[s, t] + config_param.is_Belta * sum(loads.load_curve[:, t])))
	println("\t constraints: 6) system reserves limits\t\t\t\t\t done")

	# power balance constraints
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		sum(pg₀[(1 + (s - 1) * NG):(s * NG), t]) +
		sum(winds.scenarios_curve[s, t] .* winds.p_max[:, 1] - Δpw[(1 + (s - 1) * NW):(s * NW), t]) -
		sum(loads.load_curve[:, t] - Δpd[(1 + (s - 1) * ND):(s * ND), t]).==0)

	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		sum(pg₀[(1 + (s - 1) * NG):(s * NG), t]) +
		sum(winds.scenarios_curve[s, t] * winds.p_max[:, 1] -
			Δpw[(1 + (s - 1) * NW):(s * NW), t]) - sum(loads.load_curve[:, t] - Δpd[(1 + (s - 1) * ND):(s * ND), t]).==0)
	println("\t constraints: 7) power balance constraints\t\t\t\t done")

	# ramp-up and ramp-down constraints
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		pg₀[(1 + (s - 1) * NG):(s * NG), t] -
		((t == 1) ? units.p_0[:, 1] :
		 pg₀[(1 + (s - 1) * NG):(s * NG),
			t - 1]).<=
		units.ramp_up[:, 1] .* ((t == 1) ? onoffinit[:, 1] : x[:, t - 1]) +
		units.shut_up[:, 1] .* ((t == 1) ? ones(NG, 1) : u[:, t - 1]) +
		units.p_max[:, 1] .* (ones(NG, 1) - ((t == 1) ? onoffinit[:, 1] : x[:, t - 1])))
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		((t == 1) ? units.p_0[:, 1] : pg₀[(1 + (s - 1) * NG):(s * NG), t - 1]) -
		pg₀[(1 + (s - 1) * NG):(s * NG), t].<=
		units.ramp_down[:, 1] .* x[:, t] +
		units.shut_down[:, 1] .* v[:, t] +
		units.p_max[:, 1] .* (x[:, t]))
	println("\t constraints: 8) ramp-up/ramp-down constraints\t\t\t\t done")

	# PWL constraints
	eachseqment = (units.p_max - units.p_min) / 3
	@constraint(scuc,
		[s = 1:NS, t = 1:NT, i = 1:NG],
		pg₀[i + (s - 1) * NG, t].==
		units.p_min[i, 1] * x[i, t] + sum(pgₖ[i + (s - 1) * NG, t, :]))
	@constraint(scuc,
		[s = 1:NS, t = 1:NT, i = 1:NG, k = 1:3],
		pgₖ[i + (s - 1) * NG, t, k]<=eachseqment[i, 1] * x[i, t])
	println("\t constraints: 9) piece linearization constraints\t\t\t done")

	# transmissionline power limits for basline states
	if config_param.is_NetWorkCon == 1
		for l in 1:NL
			subGsdf_units = Gsdf[l, units.locatebus]
			subGsdf_winds = Gsdf[l, winds.index]
			subGsdf_loads = Gsdf[l, loads.locatebus]
			subGsdf_psses = Gsdf[1, stroges.locatebus]
			@constraint(scuc,
				[s = 1:NS, t = 1:NT],
				sum(subGsdf_units[i] * pg₀[i + (s - 1) * NG, t] for i in 1:NG) +
				sum(subGsdf_winds[w] * (winds.scenarios_curve[s, t] * winds.p_max[w, 1] -
										Δpw[(s - 1) * NW + w, t]) for w in 1:NW) -
				sum(subGsdf_loads[d] * (loads.load_curve[d, t] - Δpd[(s - 1) * ND + d, t])
					for d in 1:ND)<=lines.p_max[l, 1])
			@constraint(scuc,
				[s = 1:NS, t = 1:NT],
				sum(subGsdf_units[i] * pg₀[i + (s - 1) * NG, t] for i in 1:NG) +
				sum(subGsdf_winds[w] * (winds.scenarios_curve[s, t] * winds.p_max[w, 1] -
										Δpw[(s - 1) * NW + w, t]) for w in 1:NW) -
				sum(subGsdf_loads[d] * (loads.load_curve[d, t] - Δpd[(s - 1) * ND + d, t])
					for d in 1:ND)>=lines.p_min[l, 1])
		end
		println("\t constraints: 10) transmissionline limits for basline\t\t\t done")
	end

	#! @𝗯𝗿𝗶𝗲𝗳: 𝘁𝗵𝗶𝘀 𝗳𝘂𝗻𝗰𝘁𝗶𝗼𝗻 𝗶𝘀 𝘂𝘀𝗲𝗱 𝘁𝗼 𝗴𝗲𝗻𝗲𝗿𝗮𝘁𝗲 𝗼𝗿𝗶𝗴𝗶𝗻𝗮𝗹 𝗻𝗲𝘁𝗹𝗼𝗮𝗱 𝗮𝗻𝗱 𝗳𝗹𝗲𝘅𝗿𝗮𝗺𝗽_𝗿𝗲𝗾𝘂𝗶𝗿𝗲𝗺𝗲𝗻𝘁
	netload_curve = zeros(NS, NT)
	for s in 1:NS
		netload_curve[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max)
	end

	flexramp_up_requirement, flexramp_dw_requirement = zeros(NS, NT), zeros(NS, NT)

	for t in 2:NT
		for s in 1:NS
			println([t s])
			flexramp_up_requirement[s, t - 1] = max(0, netload_curve[s, t] - netload_curve[s, t - 1])
			flexramp_dw_requirement[s, t - 1] = abs(min(0, netload_curve[s, t] - netload_curve[s, t - 1]))
		end
	end

	xdata = LinRange(1, NT, NT)[:, 1]
	p1 = Plots.plot(xdata, sum(loads.load_curve[d, :] for d in 1:ND), label = "load")
	p2 = Plots.plot(xdata, netload_curve[1, :], label = "netload")
	p3 = Plots.plot(xdata, flexramp_up_requirement[1, :], label = "flexramp_up_requirement")
	p4 = Plots.plot(xdata, flexramp_dw_requirement[1, :], label = "flexramp_dw_requirement")
	filepath = "D:\\IEEE_Transactions Research\\2 ieee pes\\code\\benchmark"
	p5 = Plots.plot(p1, p2, p3, p4, layout = (2, 2))
	Plots.savefig(p5, "D:/IEEE_Transactions Research/2 ieee pes/master - 1/code/benchmark/fig/re_flexramp_requirement.png")

	#! 𝗡𝗢𝗧𝗘 - 𝗳𝗹𝗲𝘅𝗿𝗮𝗺𝗽_𝗰𝗼𝗻𝘀𝘁𝗿𝗮𝗶𝗻𝘁𝘀
	κ = 1.0
	@constraint(scuc,
		[s = 1:NS, t = 1:NT],
		sum(sr⁺[(1 + (s - 1) * NG):(s * NG), t]).>=flexramp_up_requirement[s, t] * κ)

	@constraint(scuc,
		[s = 1:NS, t = 1:NT, i = 1:NG],
		sum(sr⁻[(1 + (s - 1) * NG):(s * NG), t]).>=flexramp_dw_requirement[s, t] * κ)

	# @constraint(scuc,
	# 	[s = 1:NS, t = 1:NT],
	# 	sum(units.ramp_up[g] .* x[g, t] for g in 1:NG).>=flexramp_up_requirement[s, t] * κ)

	# @constraint(scuc,
	# 	[s = 1:NS, t = 1:NT],
	# 	sum(units.ramp_down[g] .* x[g, t] for g in 1:NG).>=flexramp_dw_requirement[s, t] * κ)

	# @constraint(scuc,
	# 	[s = 1:NS, t = 1:NT],
	# 	sum(sr⁺[(1 + (s - 1) * NG):(s * NG), t]).>=sum(units.p_max[g, 1] .* u[g, t] for g in 1:NG))

	# @constraint(scuc,
	# 	[s = 1:NS, t = 1:NT],
	# 	sum(sr⁻[(1 + (s - 1) * NG):(s * NG), t]).>=sum(units.p_max[g, 1] .* v[g, t] for g in 1:NG))

	println("\n")
	println("Model has been loaded")
	println("Step-4: calculation...")
	JuMP.optimize!(scuc)

	println("callback gurobisolver\t\t\t\t\t\t\t done")
	@test JuMP.termination_status(scuc) == MOI.OPTIMAL
	@test JuMP.primal_status(scuc) == MOI.FEASIBLE_POINT
	println("#TEST: termination_status\t\t\t\t\t\t pass")

	println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
	# println("pg₀>>\n", JuMP.value.(pg₀[1:NG,1:NT]))
	# println("x>>  \n", JuMP.value.(x))
	# println("Δpd>>\n", JuMP.value.(Δpd[1:ND,1:NT]))
	# println("Δpw>>\n", JuMP.value.(Δpw[1:NW,1:NT]))
	su_cost = sum(JuMP.value.(su₀))
	sd_cost = sum(JuMP.value.(sd₀))
	pᵪ = JuMP.value.(pgₖ)
	p₀ = JuMP.value.(pg₀)
	x₀ = JuMP.value.(x)
	r⁺ = JuMP.value.(sr⁺)
	r⁻ = JuMP.value.(sr⁻)
	pᵨ = JuMP.value.(Δpd)
	pᵩ = JuMP.value.(Δpw)

	# pss_charge_state⁺ = JuMP.value.(κ⁺)
	# pss_charge_state⁻ = JuMP.value.(κ⁻)
	# pss_charge_p⁺ = JuMP.value.(pc⁺)
	# pss_charge_p⁻ = JuMP.value.(pc⁻)
	# pss_charge_q = JuMP.value.(qc)

	prod_cost = pₛ *
				c₀ *
				(sum(sum(sum(sum(pᵪ[i + (s - 1) * NG, t, :] .* eachslope[:, i] for t in 1:NT))
						 for s in 1:NS) for i in 1:NG) + sum(sum(sum(x₀[:, t] .* refcost[:, 1] for t in 1:NT)) for s in 1:NS))
	cr⁺ = pₛ *
		  c₀ *
		  sum(sum(sum(ρ⁺ * r⁺[i + (s - 1) * NG, t] for i in 1:NG) for t in 1:NT)
			  for s in 1:NS)
	cr⁻ = pₛ *
		  c₀ *
		  sum(sum(sum(ρ⁺ * r⁻[i + (s - 1) * NG, t] for i in 1:NG) for t in 1:NT)
			  for s in 1:NS)
	seq_sr⁺ = pₛ * c₀ * sum(ρ⁺ * r⁺[i, :] for i in 1:NG)
	seq_sr⁻ = pₛ * c₀ * sum(ρ⁺ * r⁻[i, :] for i in 1:NG)
	𝜟pd = pₛ * sum(sum(sum(pᵨ[(1 + (s - 1) * ND):(s * ND), t]) for t in 1:NT) for s in 1:NS)
	𝜟pw = pₛ * sum(sum(sum(pᵩ[(1 + (s - 1) * NW):(s * NW), t]) for t in 1:NT) for s in 1:NS)
	str = zeros(1, 7)
	str[1, 1] = su_cost * 10
	str[1, 2] = sd_cost * 10
	str[1, 3] = prod_cost
	str[1, 4] = cr⁺
	str[1, 5] = cr⁻
	str[1, 6] = 𝜟pd
	str[1, 7] = 𝜟pw

	open("D:\\IEEE_Transactions Research\\2 ieee pes\\master - 1\\code\\benchmark\\res\\re_flexiramp_calculation_result.txt", "w") do io
		writedlm(io, [" "])
		writedlm(io, ["su_cost" "sd_cost" "prod_cost" "cr⁺" "cr⁻" "𝜟pd" "𝜟pw"], '\t')
		writedlm(io, str, '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 1: units stutup/down states"])
		writedlm(io, JuMP.value.(x), '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 2: units dispatching power in scenario NO.1"])
		writedlm(io, JuMP.value.(pg₀[1:NG, 1:NT]), '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 3: spolied wind power"])
		writedlm(io, JuMP.value.(Δpw[1:NW, 1:NT]), '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 4: forced load curtailments"])
		writedlm(io, JuMP.value.(Δpd[1:ND, 1:NT]), '\t')
		# writedlm(io, [" "])
		# writedlm(io, ["list 5: pss charge state"])
		# writedlm(io, pss_charge_state⁺[1:NC, 1:NT], '\t')
		# writedlm(io, [" "])
		# writedlm(io, ["list 6: pss discharge state"])
		# writedlm(io, pss_charge_state⁻[1:NC, 1:NT], '\t')
		# writedlm(io, [" "])
		# writedlm(io, ["list 7: pss charge power"])
		# writedlm(io, pss_charge_p⁺[1:NC, 1:NT], '\t')
		# writedlm(io, [" "])
		# writedlm(io, ["list 8: pss discharge power"])
		# writedlm(io, pss_charge_p⁻[1:NC, 1:NT], '\t')
		# writedlm(io, [" "])
		# writedlm(io, ["list 9: pss strored energy"])
		# writedlm(io, pss_charge_q[1:NC, 1:NT], '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 10: sr⁺"])
		writedlm(io, r⁺[1:NG, 1:NT], '\t')
		writedlm(io, [" "])
		writedlm(io, ["list 11: sr⁻"])
		writedlm(io, r⁻[1:NG, 1:NT], '\t')
	end

	println("the calculation_result has been saved into | calculation_result.txt |\t done")
	println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

	# return p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, pss_charge_p⁺, pss_charge_p⁻, su_cost, sd_cost, prod_cost, cr⁺, cr⁻

	return x₀, p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cr⁺, cr⁻
end
