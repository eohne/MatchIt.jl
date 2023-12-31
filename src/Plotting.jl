"""
    balance_plot(df::MatchedIt,barmode="overlay",logscale=false)

Creates a plot with two histograms showing the balance of propensity scores before and after matching.

# Arguments:
  * df`::MatchedIt`: The output of a function call to `matchit`
  * barmode: whehter the bars of the histogram should be overlayed (`"overlay"`).
  * logscale`::Bool`: whether the yaxis should be printed on a log-scale (defaults to `false`)

# Output
A `PlotlyJS.SyncPlot`

"""
function balance_plot(df::MatchedIt,barmode="overlay",logscale=false)
    if barmode=="overlay"
        opacity = 0.8
    else
        opacity=1
    end
    fig = PlotlyJS.make_subplots(
    rows=1, cols=2,subplot_titles=["<b>Umatched</b>" "<b>Matched</b>" ])
    add_trace!(fig,
        PlotlyJS.histogram(x = df.df.Dist[isequal.(df.df[:,df.T],1)],
            name="High",opacity=opacity, marker_color="steelblue"),
        row=1,
        col=1)
    add_trace!(fig,
        PlotlyJS.histogram(x = df.df.Dist[isequal.(df.df[:,df.T],0)],
            name="Low",opacity=opacity,marker_color="orange")
        ,row=1,
        col=1)
    add_trace!(fig,
        PlotlyJS.histogram(x = df.matched.Dist[isequal.(df.matched[:,df.T],1)],
            name="High",opacity=opacity, marker_color="steelblue",showlegend=false),
        row=1,
        col=2)
    add_trace!(fig,
        PlotlyJS.histogram(x = df.matched.Dist[isequal.(df.matched[:,df.T],0)],
            name="Low",opacity=opacity,marker_color="orange",showlegend=false),
        row=1,
        col=2)
    PlotlyJS.relayout!(fig,
    barmode=barmode,
    xaxis = attr(title=attr(text="Propensity Score",standoff=10),automargin=true),
    xaxis2 = attr(title=attr(text="Propensity Score",standoff=10),automargin=true),
    plot_bgcolor=:transparent,
    paper_bgcolor = :transparent,
    yaxis = attr(title=attr(text="Frequency",standoff=10),automargin=true,gridcolor = :transparent,type=ifelse(logscale,"log",nothing)),
    yaxis2 = attr(title=attr(text="",standoff=10),automargin=true,gridcolor = :transparent,type=ifelse(logscale,"log",nothing)),
   font_family="Times New Roman",font_size=16, font_color = "black")
    return fig
end