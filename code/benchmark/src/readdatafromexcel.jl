using XLSX

function readxlssheet()
	println("Step-1: Pkgs and functions are loaded")

	df = XLSX.readxlsx("D:\\IEEE_Transactions Research\\2 ieee pes\\master - 1\\code\\benchmark\\data\\data.xlsx")
	# part-1: read frequency data
	unitsfreqparam = df["units_frequencyparam"]
	windsfreqparam = df["winds_frequencyparam"]

	Sheet1_list = string("A2", ":", "G", string(size(unitsfreqparam[:], 1)))
	Sheet2_list = string("A2", ":", "F", string(size(windsfreqparam[:], 1)))

	unitsfreqparam = convert(Array{Float64, 2}, unitsfreqparam[Sheet1_list])
	windsfreqparam = convert(Array{Float64, 2}, windsfreqparam[Sheet2_list])

	# part-2: read stroge data
	strogesystemdata = df["strogesystem_data"]
	Sheet3_list = string("A2", ":", "L", string(size(strogesystemdata[:], 1)))
	strogesystemdata = convert(Array{Float64, 2}, strogesystemdata[Sheet3_list])

	# part-3: read conventional unit, network ,and load data
	gendata = df["units_data"]
	Sheet4_list = string("A2", ":", "M", string(size(gendata[:], 1)))
	gendata = convert(Array{Float64, 2}, gendata[Sheet4_list])

	gencost = df["units_cost"]
	Sheet5_list = string("A2", ":", "H", string(size(gencost[:], 1)))
	gencost = convert(Array{Float64, 2}, gencost[Sheet5_list])

	linedata = df["branch_data"]
	Sheet6_list = string("A2", ":", "E", string(size(linedata[:], 1)))
	linedata = convert(Array{Float64, 2}, linedata[Sheet6_list])

	loadcurve = df["load_curve"]
	Sheet7_list = string("A2", ":", "B", string(size(loadcurve[:], 1)))
	loadcurve = convert(Array{Float64, 2}, loadcurve[Sheet7_list])

	loaddata = df["load_data"]
	Sheet8_list = string("A2", ":", "C", string(size(loaddata[:], 1)))
	loaddata = convert(Array{Float64, 2}, loaddata[Sheet8_list])

	return unitsfreqparam, windsfreqparam, strogesystemdata, gendata, gencost, linedata, loadcurve, loaddata
end

function forminputdata(DataGen, DataBranch, DataLoad, LoadCurve, GenCost, UnitsFreqParam, StrogeData)

	# DataGen,DataBranch,DataLoad,LoadCurve,GenCost = IEEE_RTS6()
	NB = Int64(maximum([maximum(DataBranch[:, 2]), maximum(DataBranch[:, 3])]))::Int64
	NL = Int64(size(DataBranch)[1])::Int64
	NG = Int64(size(DataGen)[1])::Int64
	ND = Int64(size(DataLoad)[1])::Int64
	NC = Int64(size(StrogeData)[1])::Int64
	NT = 24::Int64

	Gens_Index = convert(Array{Int64}, DataGen[:, 1])
	Gens_LocateBus = convert(Array{Int64}, DataGen[:, 2])
	Gens_Pmax = DataGen[:, 3] / 100
	Gens_Pmin = DataGen[:, 4] / 100
	Gens_RD = DataGen[:, 5] / 100
	Gens_RU = DataGen[:, 6] / 100
	Gens_SD = DataGen[:, 7] / 100
	Gens_SU = DataGen[:, 8] / 100
	Gens_TU = DataGen[:, 9]
	Gens_TD = DataGen[:, 10]
	Gens_x0 = DataGen[:, 11]
	# Gens_x0 = ones(NG, 1)[:, 1]
	Gens_t0 = DataGen[:, 12]
	Gens_p0 = DataGen[:, 13] / 100

	Gens_c = GenCost[:, 2]
	Gens_b = GenCost[:, 3] * 1e2
	Gens_a = GenCost[:, 4] * 1e4
	Gens_CD = GenCost[:, 5]
	Gens_CU = GenCost[:, 6]
	Gens_CU1 = GenCost[:, 7]
	Gens_Cold = GenCost[:, 8]

	Trans_index = convert(Array{Int64}, DataBranch[:, 1])
	Trans_From = convert(Array{Int64}, DataBranch[:, 2])
	Trans_To = convert(Array{Int64}, DataBranch[:, 3])
	Trans_x = DataBranch[:, 4]
	Trans_Pmax = DataBranch[:, 5] / 100
	Trans_Pmin = (-1) .* DataBranch[:, 5] / 100
	# Trans_b    = zeros(NL, 1)
	# Trans_Ratio = ones(NL, 1)

	# Loads.Curve = LoadCurve
	Loads_Index = convert(Array{Int64}, DataLoad[:, 1])
	Loads_LocateBus = convert(Array{Int64}, DataLoad[:, 2])
	Loads_Percent = DataLoad[:, 3]
	Loads_SumLoad = LoadCurve[:, 2] / 100
	Loads_PerLoad = zeros(ND, NT)
	for i in 1:NT
		Loads_PerLoad[:, i] = Loads_SumLoad[i, 1] .* Loads_Percent[:, 1]
	end

	Hg = UnitsFreqParam[:, 2]
	Dg = UnitsFreqParam[:, 3]
	Kg = UnitsFreqParam[:, 4]
	Fg = UnitsFreqParam[:, 5]
	Tg = UnitsFreqParam[:, 6]
	Rg = UnitsFreqParam[:, 7]

	Pss_index = convert(Array{Int64}, StrogeData[:, 1])
	Pss_locatebus = convert(Array{Int64}, StrogeData[:, 2])
	Pss_q_max = StrogeData[:, 3] / 100
	Pss_q_min = StrogeData[:, 4] / 100
	Pss_p⁺ = StrogeData[:, 5] / 100
	Pss_p⁻ = StrogeData[:, 6] / 100
	Pss_P₀ = StrogeData[:, 7] / 100
	Pss_γ⁺ = StrogeData[:, 8] / 100
	Pss_γ⁻ = StrogeData[:, 9] / 100
	Pss_η⁺ = StrogeData[:, 10]
	Pss_η⁻ = StrogeData[:, 11]
	Pss_δₛ = StrogeData[:, 12]

	# re-normazied data
	config_param = config(1, 1, 1, 1, 1, 3, 0.005, 0.005, 1, 1, 1, 1e5, 1e5, 50, 0.01)

	units = unit(Gens_Index, Gens_LocateBus, Gens_Pmax, Gens_Pmin, Gens_RU, Gens_RD,
		Gens_SU, Gens_SD, Gens_TU, Gens_TD, Gens_x0, Gens_t0, Gens_p0, Gens_a,
		Gens_b, Gens_c, Gens_CU, Gens_CU1, Gens_CD, Gens_Cold, Hg, Dg, Kg, Fg, Tg,
		Rg)
	# lines = transmissionline(Trans_From, Trans_To, Trans_x, Trans_b, Trans_Pmax, Trans_Pmin)

	lines = transmissionline(Trans_index, Trans_From, Trans_To, Trans_x, Trans_Pmax,
		Trans_Pmin)

	stroges = stroge(Pss_index, Pss_locatebus, Pss_q_max, Pss_q_min, Pss_p⁺, Pss_p⁻, Pss_P₀,
		Pss_γ⁺, Pss_γ⁻, Pss_η⁺, Pss_η⁻, Pss_δₛ)

	if size(Loads_PerLoad, 1) == ND
		if size(Loads_PerLoad, 2) == NT
			loads = load(Loads_Index, Loads_LocateBus, Loads_PerLoad)
		end
	end

	println("Step-2: imput data are loaded")

	return config_param, units, lines, loads, stroges, NB, NG, NL, ND, NT, NC
end
