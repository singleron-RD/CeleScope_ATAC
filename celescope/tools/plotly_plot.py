import plotly
import plotly.express as px
import plotly.figure_factory as pf


PLOTLY_CONFIG = {
    "displayModeBar": True,
    "staticPlot": False,
    "showAxisDragHandles": False,
    "modeBarButtons": [["toImage", "resetScale2d"]],
    "scrollZoom": False,
    "displaylogo": False
}


class Plotly_plot():

    def __init__(self, df):
        self._df = df
        self._fig = None

    def plotly_plot(self):
        return plotly.offline.plot(
            self._fig,
            include_plotlyjs=False,
            output_type='div',
            config=PLOTLY_CONFIG
        )

    def get_plotly_div(self):

        return self.plotly_plot()


class Insert_plot(Plotly_plot):

    def __init__(self, df):
        super().__init__(df)
        self.set_fig()
    
    def set_fig(self):

        hist_data = [self._df["isize"]]
        group_labels = ['Sample']
        colors = ['#1F77B4']

        self._fig = pf.create_distplot(hist_data, group_labels, show_hist=False, colors=colors,show_rug=False)
        self._fig.update_layout(
            title={"text": "Insert size distribution",
                   'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
            xaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            xaxis_title = "Insert Size(bp)", yaxis_title = "Density", font=dict(size=14,color="Black"),
            showlegend=False,
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF"
        )
        

class Tss_plot(Plotly_plot):

    def __init__(self, df):
        super().__init__(df)
        self.set_fig()
    
    def set_fig(self):
        self._fig = px.line(self._df, x="dists", y="agg_tss_scores", markers=True, color_discrete_sequence=['#1F77B4'])
        self._fig.update_layout(title={"text": "Aggregate TSS enrichment",
                                 'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
                          xaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                          yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
                          xaxis_title = "Distance from TSS(bp)", yaxis_title = "Aggregate TSS score", font=dict(size=14,color="Black"),
                          showlegend=False,
                          margin=dict(l=50, r=0, t=30, b=30),
                          plot_bgcolor="#FFFFFF"
        )
        
    
class Peak_plot(Plotly_plot):

    def __init__(self, df):
        super().__init__(df)
        self.set_fig()
    
    @staticmethod
    def label(df):
        if df["cell_called"] == True:
            return "Cells"
        else:
            return "Non-cells"
    
    def set_fig(self):
        self._df['cell_called'] = self._df.apply(self.label, axis=1)
        self._df.sort_values("cell_called", ascending=False, inplace=True)
        
        self._fig = px.scatter(self._df, x="total_frags", y="frac_peak", color="cell_called",color_discrete_sequence=["grey", "orange"]
                 )
        self._fig.update_layout(
            width=470, height=313,
            title={"text": "Peaks Targeting",
                   'y':0.98, 'x':0.5, 'xanchor': 'center', 'yanchor': 'top'},
            xaxis={"type": "log", "color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            xaxis_title = "Fragments per Barcode", yaxis_title = "Fraction Fragments Overlapping Peaks", font=dict(size=12,color="Black"),
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF",
            legend=dict(title_text='')
        )
        self._fig.update_traces(marker_size=3)