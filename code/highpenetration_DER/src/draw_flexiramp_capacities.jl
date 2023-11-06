function draw_flexiramp(flexiramp_up_requirement_1, flexiramp_up_capacity_1, flexiramp_dw_requirement_1, flexiramp_dw_capacity_1, selected_scenario, NT)
    default(legend_font_halign = :left)

    ctg = repeat([L"\textrm{online\,\, flexiramp\,\, requirement}\,\,\uparrow", L"\textrm{online\,\, flexiramp\,\, capacity}\,\,\uparrow"], inner = NT)
    nam = repeat(collect(1:1:NT), outer = 2)
    p4 = groupedbar(nam, hcat(flexiramp_up_requirement_1[selected_scenario, :], flexiramp_up_capacity_1[selected_scenario, :]),
        group = ctg,
        size = (300, 280),
        ylims = (-1, 4.5),
        grid = false,
        legendfontsize = 8,
        tickfontfamily = "Palatino Bold",
        # tickfontfamily = "Palatino Bold",
        xlabel = L"t / h",
        # fillcolor= [:orange, :blue],
        ylabel = L"\textrm{Power / p.u.}",
        bar_width = 0.55,
        # color = :cyclic_grey_15_85_c0_n256,
        lw = 0,
        framestyle = :box)

    ctg = repeat([L"\textrm{online\,\, flexiramp\,\, requirement}\,\,\downarrow", L"\textrm{online\,\, flexiramp\,\, capacity}\,\,\downarrow"], inner = NT)
    p4 = groupedbar!(nam, -hcat(flexiramp_dw_requirement_1[selected_scenario, :], flexiramp_dw_capacity_1[selected_scenario, :]),
        group = ctg,
        size = (300, 280),
        # ylims = (0, 3),
        grid = false,
        # fillcolor= [:orange, :blue],
        # legendfontfamily = "Palatino Bold",
        # tickfontfamily = "Palatino Bold",
        xlabel = L"t / h",
        ylabel = L"\textrm{Power / p.u.}",
        bar_width = 0.55,
        lw = 0,
        framestyle = :box)
    return p4
end
