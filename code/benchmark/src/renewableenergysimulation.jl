# windsimulation
using Distributions

struct wind
    index::Vector{Int64}
    locatebus::Vector{Int64}
    p_max::Vector{Float64}
    scenarios_prob::Float64
    scenarios_nums::Int64
    scenarios_curve::Array{Float64}

    # frequency constrol process of renewable energy
    Fcmode::Vector{Float64}
    # model-1(drop constrol):Kw/(Rw * (1 + sTw))
    Kw::Vector{Float64}
    Rw::Vector{Float64}
    # model-2(Virtual inertia mechine): (sMw + Dw)/(1 + sTw)
    Mw::Vector{Float64}
    Dw::Vector{Float64}
    Tw::Vector{Float64}
    # wind(Fcmode,Kw,Rw,Mw,Dw) = new(Fcmode,Kw,Rw,Mw,Dw)
    function checkvaildity(Fcmode)
        if Fcmode == 1
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-1 and related params")
            end
        end
        return if Fcmode == 2
            if Kw != 0 || Rw != 0
                println("dismathch for Fcmode-2 and related params")
            end
        end
    end
    function wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve,
                  Fcmode, Kw, Rw, Mw, Dw, Tw)
        return new(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve,
                   Fcmode, Kw, Rw, Mw, Dw, Tw)
    end
end

# predefine value
index = [1; 2]
locatebus = [1; 1]
NW = length(index)
scenarios_nums = 5
scenarios_prob = 1 / scenarios_nums
# scenarios_curve = zeros(scenarios_nums, NT)
NT = 24
# assum the capacity of each wind is same and reforced as 0.5 p.u.
cap = [0.5] * 1.0
p_max = cap .* ones(NW, 1)
p_max = p_max[:, 1]
scenarios_curve = zeros(scenarios_nums, NT)
scenarios_curvebase = [
    0.440724927203680 0.420965256587272 0.449034794022911 0.454128108336623 0.436483077739172 0.477450522402300
    0.443871634609799 0.374756446192485 0.448192193924943 0.431190577826877 0.428867647037057 0.445673091565042
    0.433764408789611 0.421900481861469 0.429104412188035 0.463277796146724 0.426579282372516 0.448189506134410
    0.429353980231385 0.434861266141317 0.437494540514197 0.456877055120346 0.425139803090161 0.425629623577982
]
scenarios_curvebase = reshape(scenarios_curvebase, 1, NT)

function genscenario(WindsFreqParam)
    # srand(123)
    #     sample_sets = rand(Weibull(), scenarios_nums * NT) * 0.01
    #     scenarios_error = reshape(sample_sets, scenarios_nums, NT)
    #     for i = 1 : scenarios_nums
    #         for j in 1 : NT
    #             sample_temp = rand()
    #             if sample_temp > 0.5
    #                 scenarios_curve[i, j] = scenarios_curvebase[1, j] + scenarios_error[i, j]
    #             else
    #                 scenarios_curve[i, j] = scenarios_curvebase[1, j] - scenarios_error[i, j]
    #             end
    #         end
    #     end
    scenarios_curve = [
        0.440724927203680 0.420965256587272 0.449034794022911 0.454128108336623 0.436483077739172 0.477450522402300 0.443871634609799 0.374756446192485 0.448192193924943 0.431190577826877 0.428867647037057 0.445673091565042 0.433764408789611 0.421900481861469 0.429104412188035 0.463277796146724 0.426579282372516 0.448189506134410 0.429353980231385 0.434861266141317 0.437494540514197 0.456877055120346 0.425139803090161 0.425629623577982
        0.438145251438362 0.451595831499290 0.434476599311993 0.419306858427854 0.439299123016117 0.402675152643531 0.436348294887821 0.447513027575036 0.445276832579360 0.408448500875771 0.476106019486472 0.451932867123187 0.446968204950444 0.457706023689642 0.454429491703142 0.432489551344388 0.460269791720502 0.417994780067730 0.404420416693225 0.443013967794901 0.407382847053778 0.430503777173583 0.455183618944849 0.443789093804304
        0.441672518788126 0.461922845597782 0.425338820890952 0.420366090471607 0.411612893296905 0.435840069094316 0.443930499695973 0.457511112047526 0.450817300160177 0.396413160907573 0.441068179613219 0.432401166117165 0.420639320150678 0.443835529493502 0.433537192471826 0.427399090347307 0.417573417186437 0.422905158624658 0.467119379846108 0.500219495833784 0.432716758754643 0.422622895486611 0.452734491278974 0.425638923917095
        0.404454468849607 0.427898443023513 0.456678201304931 0.466227716201764 0.458275535897303 0.447722201714402 0.430408791524416 0.457075404497126 0.422560643262941 0.479016292670867 0.440735317418680 0.426724520859767 0.438876736399005 0.427232619437777 0.431855284702541 0.436026678964276 0.463231975065251 0.449494208561186 0.440539555865053 0.430643201150598 0.463252275085078 0.426482114000725 0.450319628567473 0.447144352980500
        0.475002833720225 0.437617623292143 0.434471584469213 0.439971226562153 0.454329370050504 0.436312054145452 0.445440779281991 0.463144009687827 0.433153030072578 0.484931467718911 0.413222836444571 0.443268354334839 0.459751329710262 0.449325345517610 0.451073618934457 0.440806883197305 0.432345533655294 0.461416346612016 0.458566667364229 0.391262069079400 0.459153578592304 0.463514158218735 0.416622458118542 0.457798005720119
    ]

    FCmode = WindsFreqParam[:, 1]
    KW = WindsFreqParam[:, 2]
    RW = WindsFreqParam[:, 3]
    MW = WindsFreqParam[:, 4]
    DW = WindsFreqParam[:, 5]
    TW = WindsFreqParam[:, 6]

    winds = wind(index, locatebus, p_max, scenarios_prob, scenarios_nums, scenarios_curve,
                 FCmode, KW, RW, MW, DW, TW)

    return winds, NW
end
