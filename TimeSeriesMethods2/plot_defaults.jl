using PlotThemes
using PlotUtils

mints_palette = [
    colorant"#3cd184", # mint green
    colorant"#f97171", # dark coral
    colorant"#1e81b0", # dark blue
    colorant"#66beb2", # dark blue-green
    colorant"#f99192", # light coral
    colorant"#8ad6cc", # middle blue-green
    colorant"#3d6647", # dark green
    #        colorant"#8FDDDF", # middle blue
]



function add_mints_theme()
    _mints = PlotThemes.PlotTheme(Dict([
        # :background => :white,
        :framestyle => :box,
        :tickdirection => :out,
        :linewidth => 2,
        :markerstrokewidth => 0,
        :colorgradient => :inferno, #:vik,  # or :turbo for Dr. Lary's preference
        :titlefontsize => 15,
        :labelfontsize => 13,
        :tickfontsize => 11,
        :colorbar_titlefontsize => 13,
        :palette => mints_palette,
        :grid => true,
        :minorgrid => true,
        :gridwidth => 2,
        :minorgridwidth => 1,
        :gridalpha => 0.3,
        :minorgridalpha => 0.25,
        #:minorticks => 5,
        :legend => :outertopright,
    ]))

    PlotThemes.add_theme(:mints, _mints)
end


# myfont = Plots.font("sans-serif", :black)

# default(
#     linewidth=2,
#     grid=true,
#     tickdirection=:out,
#     minorgrid=true,
#     gridwidth=2,
#     minorgridwidth=1,
#     gridalpha=0.3,
#     minorgridalpha=0.25,
#     margin=5*Plots.mm,
#     topmargin=2*Plots.mm,
#     framestyle=:box,
#     titlefont=myfont,
#     guidefont=myfont,
#     legendfont=myfont,
#     tickfont=myfont,
#     titlefontsize=16,
#     labelfontsize=13,
#     tickfontsize=11,
#     colorbar_titlefontsize=13,
#     overwrite_figure=false,
# )

