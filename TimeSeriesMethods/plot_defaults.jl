myfont = Plots.font("sans-serif", :black)

default(
    linewidth=2,
    grid=true,
    tickdirection=:out,
    minorgrid=true,
    gridwidth=2,
    minorgridwidth=1,
    gridalpha=0.3,
    minorgridalpha=0.25,
    margin=5*Plots.mm,
    topmargin=2*Plots.mm,
    framestyle=:box,
    titlefont=myfont,
    guidefont=myfont,
    legendfont=myfont,
    tickfont=myfont,
    titlefontsize=16,
    labelfontsize=13,
    tickfontsize=11,
    colorbar_titlefontsize=13,
    overwrite_figure=false,
)

