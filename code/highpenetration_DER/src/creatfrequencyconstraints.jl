using MultivariateStats

function creatfrequencyfittingfunction(units, winds, NG, NW)
	Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, sampleStatues = montecalrosimulation(units, winds, NG, NW)

	# x: |---H---|---D---|---F---|---K---|
	# y: f_nadir
	y = Set_f_nadir
	x = zeros(size(set_H, 1), 5)
	for i in 1:size(set_H, 1)
		x[i, 1:4] = [set_H[i, 1], set_Dg[i, 1], set_Fg[i, 1], set_Kg[i, 1]]
		if isnan(set_H[i, 1]) == true
			x[i, 2], x[i, 3], x[i, 4], x[i, 5] = NaN, NaN, NaN, NaN
		else
			x[i, 5] = set_δp[i, 1]
		end
	end

	# clearing the obtained data
	# x₀ = filter(!isnan, x)
	# y₀ = filter(!isnan, y)
	x₀ = hcat((filter(!isnan, x[:, i]) for i in 1:size(x, 2))...)
	y₀ = hcat((filter(!isnan, y[:, i]) for i in 1:size(y, 2))...)

	# delate the duplicate elements
	# x₁ = Rmduplicateelements(x₀, 1)
	# y₁ = Rmduplicateelements(y₀, 1)

	# println(x), println(y)
	sol = llsq(x₀, y₀)
	A, b = sol[1:(end - 1), :], sol[end, :]

	return A, b
end

function montecalrosimulation(units, winds, NG, NW)
	# create unitsstatues
	sampleNumber = 200
	unitsamplestatues = rand(0:1, NG, sampleNumber)
	Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, set_Rg = zeros(sampleNumber, 1),
	zeros(sampleNumber, 1),
	zeros(sampleNumber, 1),
	zeros(sampleNumber, 1),
	zeros(sampleNumber, 1),
	zeros(sampleNumber, 1),
	zeros(sampleNumber, 1)
	for i in 1:sampleNumber
		Sampling_Statue = unitsamplestatues[:, i]
		f_nadir, t_nadir, H, δp, Kg, Fg, Rg, Dg = creatingfrequencyresponsesamplingdata(units, winds, NW, NG, Sampling_Statue)
		Set_f_nadir[i, 1], set_H[i, 1], set_δp[i, 1], set_Dg[i, 1], set_Fg[i, 1], set_Kg[i, 1], set_Rg[i, 1] = f_nadir, H, δp, Dg, Fg, Kg, Rg
	end

	return Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, unitsamplestatues
end

function creatingfrequencyresponsesamplingdata(units, winds, NW, NG, Sampling_Statue)

	# normalized winds parameters through COI
	vsmFC_number = sum(winds.Fcmode[:, 1])
	doopFC_number = length(winds.Fcmode[:, 1]) - vsmFC_number
	adjustablewindsVSCpower = winds.Fcmode .* winds.p_max
	# current_Kw    = sum(winds.Kw ./ winds.Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max)) # Kw
	inverse_winds_Rw = zeros(NW, 1)
	for i in 1:NW
		if winds.Fcmode[i, 1] == 0
			inverse_winds_Rw[i, 1] = 1 / winds.Rw[i, 1]
		end
	end
	current_Kw = sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max))
	current_Dw = sum(winds.Dw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Dw
	current_Mw = sum(winds.Mw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Mw
	current_Hw = current_Mw / 2
	current_Rw = 1.0

	# units parameters
	adjustabletheramlpower = units.p_max .* Sampling_Statue
	current_Kg = sum(units.Kg ./ units.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # Kg
	current_Tg = sum(units.Tg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # T
	current_Fg_div_Rg = sum(units.Kg .* units.Fg ./ units.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
	current_Fg = current_Fg_div_Rg
	current_Rg = 1
	current_Dg = sum(units.Dg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
	current_Hg = sum(units.Hg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
	current_Mg = current_Hg * 2

	# total powers
	#  powers for intia frequency response
	localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
	#  powers for primary and second frequency responses
	sumapparentpower = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
	# potential prob. # δp
	p_step = maximum(units.p_max)

	# sumD and sumH
	current_sumD = (sum(units.Dg .* adjustabletheramlpower) +
					sum(winds.Dw .* winds.p_max .* winds.Fcmode + winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode))) / sumapparentpower
	current_sumH = (sum(current_Mg .* adjustabletheramlpower) +
					sum(current_Mw .* adjustablewindsVSCpower)) / (sum(adjustabletheramlpower) + sum(adjustablewindsVSCpower)) / 2

	# f_db = 0.2 # param-1 of deadarea frequencychange for units
	# t_db = -1 * current_M / current_D * log(1 - f_db * current_D / (current_M * p_step))

	# simplify parameters
	# part-1 units parameters
	Mg, Hg, Dg, Tg, Rg, Fg, Kg = current_Mg, current_Hg, current_Dg, current_Tg, current_Rg, current_Fg, current_Kg
	# part-2 winds parameters
	Rw, Dw, Mw, Hw, Tw = current_Rw, current_Dw, current_Mw, current_Hw, current_Tg
	# part-3 power unbalances
	δp = p_step
	D, H, T = current_sumD, current_sumH, Tg

	# a1, a2, a3, a4, a₅
	a₁ = 2 * H * T * Rg * Rw
	a₂ = 0.5 * Rg * Rw * (2 * H + T * D) + (Mw * Rg * Rw + Fg * Kg * Rw * T)
	a₃ = Rw * Rg + Kg * Rw + Dw * Rg * Rw
	a₄ = Rg * Rw / a₁
	a₅ = a₄ * δp

	β₁ = a₁ / a₃
	γ₁ = sqrt(a₂^2 - a₁ * a₃)
	t_nadir = a₁ * sqrt(a₂^2 - a₁ * a₃) * acosh(abs((a₁ - a₂ * T) / (sqrt(a₁) * sqrt(abs(a₁ - 2 * a₂ * T + a₃ * (T^2)))))) / γ₁
	ℓ₁ = γ₁ / a₁
	Δf = β₁ * (exp(-a₂ / a₁) * cosh(ℓ₁ * t_nadir / 10) - 1 + (a₂ - a₃ * T) / ℓ₁ * sinh(ℓ₁ * t_nadir / 10)) * (-1) / 50
	f_nadir = 50 - abs(Δf)

	return f_nadir, t_nadir, H, δp, Kg, Fg, Rg, Dg
end

# delate the duplicate elements in arrays by filitering index
function Rmduplicateelements(str, index)
	str_index = unique(str[:, index])
	str_1 = zeros(size(str_index, 1), size(str, 2))
	for i in 1:size(str_index, 1)
		j = findfirst(isequal(str_index[i, 1]), str[:, index])
		str_1[i, :] = str[j, :]
	end

	return str_1
end

# function creatfrequencyfittingfunction(units, winds, NG, NW)

#     Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, sampleStatues = montecalrosimulation(units, winds, NG, NW)

#     # x: |---H---|---D---|---F---|---K---|
#     # y: f_nadir
#     y = Set_f_nadir
#     x = zeros(size(set_H, 1), 5)
#     for i in 1:size(set_H, 1)
#         x[i,1:4] = [set_H[i,1], set_Dg[i,1], set_Fg[i,1], set_Kg[i,1]]
#         if isnan(set_H[i,1]) == true
#             x[i,2], x[i,3], x[i,4], x[i,5] = NaN, NaN, NaN, NaN
#         else
#             x[i,5] = set_δp[i,1]
#         end
#     end

#     # clearing the obtained data
#     # x₀ = filter(!isnan, x)
#     # y₀ = filter(!isnan, y)
#     x₀ = hcat((filter(!isnan, x[:, i]) for i in 1:size(x, 2))...)
#     y₀ = hcat((filter(!isnan, y[:, i]) for i in 1:size(y, 2))...)

#     # delate the duplicate elements
#     # x₁ = Rmduplicateelements(x₀, 1)
#     # y₁ = Rmduplicateelements(y₀, 1)

#     # println(x), println(y)
#     sol = llsq(x₀, y₀)
#     A, b = sol[1:end - 1,:], sol[end,:]

#     return A, b

# end

# function montecalrosimulation(units, winds, NG, NW)
#     # create unitsstatues
#     sampleNumber = 200
#     unitsamplestatues = rand(0:1, NG, sampleNumber)
#     Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg = zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1), zeros(sampleNumber, 1)
#     for i in 1:sampleNumber
#         Sampling_Statue = unitsamplestatues[:,i]
#         f_nadir, t_nadir, H, δp, Fg, Kg, Fg = creatingfrequencyresponsesamplingdata(units, winds, NW, NG, Sampling_Statue)
#         Set_f_nadir[i,1], set_H[i,1], set_δp[i,1], set_Dg[i,1], set_Fg[i,1], set_Kg[i,1] = f_nadir, t_nadir, H, δp, Fg, Kg, Fg
#     end

#     return Set_f_nadir, set_H, set_δp, set_Dg, set_Fg, set_Kg, unitsamplestatues

# end

# function creatingfrequencyresponsesamplingdata(units, winds, NW, NG, Sampling_Statue)

#     # normalized winds parameters through COI
#     vsmFC_number  = sum(winds.Fcmode[:, 1])
#     doopFC_number = length(winds.Fcmode[:, 1]) - vsmFC_number
#     adjustablewindsVSCpower = winds.Fcmode .* winds.p_max
#     # current_Kw    = sum(winds.Kw ./ winds.Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max)) # Kw
#     inverse_winds_Rw = zeros(NW, 1)
#     for i in 1:NW
#         if winds.Fcmode[i,1] == 0
#             inverse_winds_Rw[i,1] = 1 / winds.Rw[i,1]
#         end
#     end
#     current_Kw    = sum(winds.Kw .* inverse_winds_Rw .* (ones(NW, 1) - winds.Fcmode) .* winds.p_max) / sum(((ones(NW, 1) - winds.Fcmode) .* winds.p_max))
#     current_Dw    = sum(winds.Dw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Dw
#     current_Mw    = sum(winds.Mw .* adjustablewindsVSCpower) / sum(adjustablewindsVSCpower) # Mw
#     current_Hw    = current_Mw / 2
#     current_Rw    = 1.0

#     # units parameters
#     adjustabletheramlpower = units.p_max .* Sampling_Statue
#     current_Kg        = sum(units.Kg ./ units.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # Kg
#     current_Tg        = sum(units.Tg .* adjustabletheramlpower) / sum(adjustabletheramlpower) # T
#     current_Fg_div_Rg = sum(units.Kg .* units.Fg ./ units.Rg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
#     current_Fg        = current_Fg_div_Rg
#     current_Rg        = 1
#     current_Dg        = sum(units.Dg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
#     current_Hg        = sum(units.Hg .* adjustabletheramlpower) / sum(adjustabletheramlpower)
#     current_Mg        = current_Hg * 2

#     # total powers
#     #  powers for intia frequency response
#     localapparentpower = (sum(units.p_max[:, 1]) + sum(winds.p_max .* winds.Fcmode))
#     #  powers for primary and second frequency responses
#     sumapparentpower   = (localapparentpower - sum(winds.p_max .* winds.Fcmode) + sum(winds.p_max))
#     # potential prob. # δp
#     p_step             = maximum(units.p_max)

#     # sumD and sumH
#     current_sumD = (
#         sum(units.Dg .* adjustabletheramlpower) +
#         sum(winds.Dw .* winds.p_max .* winds.Fcmode + winds.Kw .* winds.p_max .* (ones(NW, 1) - winds.Fcmode))
#         ) / sumapparentpower
#     current_sumH = (
#         sum(current_Mg .* adjustabletheramlpower) +
#         sum(current_Mw .* adjustablewindsVSCpower)
#         ) / (sum(adjustabletheramlpower) + sum(adjustablewindsVSCpower)) / 2

#     # f_db = 0.2 # param-1 of deadarea frequencychange for units
#     # t_db = -1 * current_M / current_D * log(1 - f_db * current_D / (current_M * p_step))

#     # simplify parameters
#     # part-1 units parameters
#     Mg, Hg, Dg, Tg, Rg, Fg, Kg = current_Mg, current_Hg, current_Dg, current_Tg, current_Rg, current_Fg, current_Kg
#     # part-2 winds parameters
#     Rw, Dw, Mw, Hw, Tw = current_Rw, current_Dw, current_Mw, current_Hw, current_Tg
#     # part-3 power unbalances
#     δp = p_step
#     D, H, T = current_sumD, current_sumH, Tg

#     # a1, a2, a3, a4, a₅
#     a₁ = 2 * H * T * Rg * Rw
#     a₂ = 0.5 * Rg * Rw * (2 * H + T * D) + (Mw * Rg * Rw + Fg * Kg * Rw * T)
#     a₃ = Rw * Rg + Kg * Rw + Dw * Rg * Rw
#     a₄ = Rg * Rw / a₁
#     a₅ = a₄ * δp

#     β₁ = a₁ / a₃
#     γ₁ = sqrt(a₂^2 - a₁ * a₃)
#     t_nadir = a₁ * sqrt(a₂^2 - a₁ * a₃) * acosh(abs((a₁ - a₂ * T) / (sqrt(a₁) * sqrt(abs(a₁ - 2 * a₂ * T + a₃ * (T^2)))))) / γ₁
#     ℓ₁ = γ₁ / a₁
#     f_nadir = β₁ * (exp(-a₂ / a₁) * cosh(ℓ₁ * t_nadir) - 1 + (a₂ - a₃ * T) / ℓ₁ * sinh(ℓ₁ * t_nadir)) * (-1)
#     f_nadir = abs(f_nadir)

#     return f_nadir, t_nadir, H, δp, Fg, Kg, Fg

# end

# # delate the duplicate elements in arrays by filitering index
# function Rmduplicateelements(str, index)
#     str_index = unique(str[:, index])
#     str_1  = zeros(size(str_index, 1), size(str, 2))
#     for i in 1:size(str_index, 1)
#         j = findfirst(isequal(str_index[i,1]), str[:,index])
#         str_1[i,:] = str[j,:]
#     end

#     return str_1

# end
