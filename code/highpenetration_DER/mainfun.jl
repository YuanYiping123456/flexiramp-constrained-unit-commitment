#NOTE: Code by Yiping Yuan

using Revise, JuMP, Gurobi, Test, DelimitedFiles, Plots, PlotlyJS, LaTeXStrings
using StatsPlots, RecipesPipeline
include("src/formatteddata.jl")
include("src/renewableenergysimulation.jl")
include("src/showboundrycase.jl")
include("src/readdatafromexcel.jl")
include("src/uccommitmentmodel.jl")
include("src/casesploting.jl")
include("src/creatfrequencyconstraints.jl")
include("src/flexramp_uccommitmentmodel.jl")
include("src/re_flexramp_uccommitmentmodel.jl")
include("src/draw_UpandDownforward_reserve.jl")
include("src/calculate_netload_and_flexiramp.jl")
include("src/calculate_flexiramp.jl")
include("src/draw_flexiramp_capacities.jl")

UnitsFreqParam, WindsFreqParam, StrogeData, DataGen, GenCost, DataBranch, LoadCurve, DataLoad = readxlssheet()
config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC = forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

if config_param.is_WindIntegration == 1
	winds, NW = genscenario(WindsFreqParam)
	# for i in 2:NT
	#     oplot(x,winds.scenarios_curve[i,:])
	# end
	Plots.plot(winds.scenarios_curve')
end
penetration_rate = sum(winds.p_max) / (sum(winds.p_max) + sum(units.p_max))
penetration_rate = sum(winds.p_max) / maximum(sum(loads.load_curve, dims = 1)[1, :])
# boundrycondition(NB::Int64,NL::Int64,NG::Int64,NT::Int64,ND::Int64,units::unit,loads::load,lines::transmissionline,winds::wind,stroges::stroge)
#! construct uc and redesigned ucs

ben_x₀, ben_p₀, ben_pᵨ, ben_pᵩ, ben_seq_sr⁺, ben_seq_sr⁻, ben_su_cost, ben_sd_cost, ben_prod_cost, ben_cost_sr⁺, ben_cost_sr⁻ = scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)
x₀, p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻ = flexramp_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)
re_x₀, re_p₀, re_pᵨ, re_pᵩ, re_seq_sr⁺, re_seq_sr⁻, re_su_cost, re_sd_cost, re_prod_cost, re_cost_sr⁺, re_cost_sr⁻ = re_flexramp_scucmodel(NT, NB, NG, ND, NC, units, loads, winds, lines, config_param)

# plotcasestudies(p₀, pᵨ, pᵩ, seq_sr⁺, seq_sr⁻, su_cost, sd_cost, prod_cost, cost_sr⁺, cost_sr⁻, NT, NG, ND, NW, NC)
#? ---figure
p1 = draw_UpandDownforward_reserve(ben_seq_sr⁺, ben_seq_sr⁻, seq_sr⁺, seq_sr⁻, re_seq_sr⁺, re_seq_sr⁻)
# netload_curve_1, netload_curve_2, netload_curve_3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
# tem1, tem2, tem3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
# str1, str2, str3 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
# for s in 1:NS
# 	t1 = (s - 1) * NW + 1
# 	t2 = s * NW
# 	tem1[s, :] .= sum(ben_pᵩ[t1:t2, :], dims = 1)[1, :]
# 	tem2[s, :] .= sum(pᵩ[t1:t2, :], dims = 1)[1, :]
# 	tem3[s, :] .= sum(re_pᵩ[t1:t2, :], dims = 1)[1, :]
# 	str1[s, :] .= sum(pᵨ[t1:t2, :], dims = 1)[1, :]
# 	str2[s, :] .= sum(re_pᵨ[t1:t2, :], dims = 1)[1, :]
# 	str3[s, :] .= sum(ben_pᵨ[t1:t2, :], dims = 1)[1, :]
# end

# for s in 1:NS
# 	netload_curve_1[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem1[s, :] - str1[s, :]
# 	netload_curve_2[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem2[s, :] - str2[s, :]
# 	netload_curve_3[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[s, :] .* sum(winds.p_max) + tem3[s, :] - str3[s, :]
# end

# netload_curve_4, netload_curve_5, netload_curve_6 = zeros(NS, NT), zeros(NS, NT), zeros(NS, NT)
# for s in 1:NS
# 	netload_curve_4[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem1[s, :] - str1[s, :]
# 	netload_curve_5[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[1, :] .* sum(winds.p_max) + tem2[s, :] - str2[s, :]
# 	netload_curve_6[s, :] .= sum(loads.load_curve[d, :] for d in 1:ND) - winds.scenarios_curve[s, :] .* sum(winds.p_max) + tem3[s, :] - str3[s, :]
# end

# Plots.plot(netload_curve_1[1, :])
# Plots.plot!(netload_curve_2[1, :])
# Plots.plot!(netload_curve_3[1, :])
# Plots.plot!(netload_curve_4[1, :], lw = 4)
# Plots.plot!(netload_curve_5[1, :], lw = 4)
# Plots.plot!(netload_curve_6[1, :], lw = 4)

#! @𝗯𝗿𝗶𝗲𝗳: 𝘁𝗵𝗶𝘀 𝗳𝘂𝗻𝗰𝘁𝗶𝗼𝗻 𝗶𝘀 𝘂𝘀𝗲𝗱 𝘁𝗼 𝗴𝗲𝗻𝗲𝗿𝗮𝘁𝗲 𝗼𝗿𝗶𝗴𝗶𝗻𝗮𝗹 𝗻𝗲𝘁𝗹𝗼𝗮𝗱 𝗮𝗻𝗱 𝗳𝗹𝗲𝘅𝗿𝗮𝗺𝗽_𝗿𝗲𝗾𝘂𝗶𝗿𝗲𝗺𝗲𝗻𝘁
NS = 5
netload_curve_1, netload_curve_2, netload_curve_3 = calculate_netload_curve(ben_pᵩ, pᵩ, re_pᵩ, ben_pᵨ, pᵨ, re_pᵨ, loads, winds)
flexiramp_up_requirement_1, flexiramp_dw_requirement_1 = calculate_UpandDw_flexramp_requirement(netload_curve_1, NT, NS)
flexiramp_up_requirement_2, flexiramp_dw_requirement_2 = calculate_UpandDw_flexramp_requirement(netload_curve_2, NT, NS)
flexiramp_up_requirement_3, flexiramp_dw_requirement_3 = calculate_UpandDw_flexramp_requirement(netload_curve_3, NT, NS)

flexiramp_up_capacity_1, flexiramp_dw_capacity_1 = calculate_flexiramp_capacity(units, ben_x₀, ben_p₀, NS, NG, NT)
flexiramp_up_capacity_2, flexiramp_dw_capacity_2 = calculate_flexiramp_capacity(units, x₀, p₀, NS, NG, NT)
flexiramp_up_capacity_3, flexiramp_dw_capacity_3 = calculate_flexiramp_capacity(units, re_x₀, re_p₀, NS, NG, NT)

# Plots.plot(collect(1:1:NT), zeros(1,NT)[1,:], fillrange = flexiramp_up_capacity_1[s,:], fillalpha = 0.35, c = 1, label = "pro-defined upforward flexiramp")
# Plots.plot!(collect(1:1:NT), zeros(1,NT)[1,:], fillrange = flexiramp_dw_capacity_1[s,:], fillalpha = 0.35, c = 1, label = "pro-defined upforward flexiramp")
# Plots.plot!(flexiramp_up_requirement_1[s, :])
# Plots.plot(flexiramp_up_capacity_2[s,:])
# Plots.plot!(flexiramp_up_requirement_2[s, :])
# Plots.plot(flexiramp_up_capacity_3[s,:])
# Plots.plot!(flexiramp_up_requirement_3[s, :])

netload_curve_4, netload_curve_5, netload_curve_6 = calculate_ini_netload_curve(ben_pᵩ, pᵩ, re_pᵩ, ben_pᵨ, pᵨ, re_pᵨ, loads, winds)
flexiramp_up_requirement_4, flexiramp_dw_requirement_4 = calculate_UpandDw_flexramp_requirement(netload_curve_4, NT, NS)
flexiramp_up_requirement_5, flexiramp_dw_requirement_5 = calculate_UpandDw_flexramp_requirement(netload_curve_5, NT, NS)
flexiramp_up_requirement_6, flexiramp_dw_requirement_6 = calculate_UpandDw_flexramp_requirement(netload_curve_6, NT, NS)
flexiramp_up_capacity_4, flexiramp_dw_capacity_4 = calculate_flexiramp_capacity(units, ben_x₀, ben_p₀, NS, NG, NT)
flexiramp_up_capacity_5, flexiramp_dw_capacity_5 = calculate_flexiramp_capacity(units, x₀, p₀, NS, NG, NT)
flexiramp_up_capacity_6, flexiramp_dw_capacity_6 = calculate_flexiramp_capacity(units, re_x₀, re_p₀, NS, NG, NT)
selected_scenario = 2

selected_scenario = 2
p2 = draw_flexiramp(flexiramp_up_requirement_1, flexiramp_up_capacity_1, flexiramp_dw_requirement_1, flexiramp_dw_capacity_1, selected_scenario, NT)
p3 = draw_flexiramp(flexiramp_up_requirement_2, flexiramp_up_capacity_2, flexiramp_dw_requirement_2, flexiramp_dw_capacity_2, selected_scenario, NT)
p3 = Plots.plot!(collect(1:1:NT),
	flexiramp_up_requirement_5[selected_scenario, :],
	lw = 0.75,
	lc = :red,
	markerstrokewidth = 0.75,
	markerstrokecolor = :blue,
	markershape = :diamond,
	markersize = 2,
	label = L"\textrm{predefined\,\, flexiramp\,\, requirement}\,\,\uparrow")
p3 = Plots.plot!(collect(1:1:NT),
	-flexiramp_dw_requirement_5[selected_scenario, :],
	lw = 0.75,
	lc = :black,
	markerstrokewidth = 0.75,
	markerstrokecolor = :blue,
	markershape = :diamond,
	markersize = 2,
	label = L"\textrm{predefined\,\, flexiramp\,\, requirement}\,\,\downarrow")
p4 = draw_flexiramp(flexiramp_up_requirement_3, flexiramp_up_capacity_3, flexiramp_dw_requirement_3, flexiramp_dw_capacity_3, selected_scenario, NT)

# p4 = Plots.plot!(collect(1:1:NT), flexramp_up_requirement[1, :])
# p4 = Plots.plot!(collect(1:1:NT), seq_sr⁺)

p5 = Plots.plot(p2, p3, p4, size = (300, 840), layout = (3, 1))
filepath = pwd()
Plots.savefig(p5, filepath * "\\highpenetration_DER\\fig\\flexiramp_compare.svg")
Plots.savefig(p1, filepath * "\\highpenetration_DER\\fig\\flexiramp_compare.pdf")

p6 = Plots.plot(p3, p4, size = (600, 280), layout = (1, 2))
Plots.savefig(p6, filepath * "\\highpenetration_DER\\fig\\flexiramp_compare_2.svg")
Plots.savefig(p6, filepath * "\\highpenetration_DER\\fig\\flexiramp_compare_2.pdf")
# p5 = draw_flexiramp(flexiramp_up_requirement_4, flexiramp_up_capacity_4, flexiramp_dw_requirement_4, flexiramp_dw_capacity_4, selected_scenario, NT)
# p6 = draw_flexiramp(flexiramp_up_requirement_5, flexiramp_up_capacity_5, flexiramp_dw_requirement_5, flexiramp_dw_capacity_5, selected_scenario, NT)
# p7 = draw_flexiramp(flexiramp_up_requirement_6, flexiramp_up_capacity_6, flexiramp_dw_requirement_6, flexiramp_dw_capacity_6, selected_scenario, NT)
