import plotly
import numpy as np
import plotly.express as px
import plotly.figure_factory as pf
import plotly.graph_objects as go
from collections import defaultdict
from celescope.tools import utils


PLOTLY_CONFIG = {
    "displayModeBar": True,
    "staticPlot": False,
    "showAxisDragHandles": False,
    "modeBarButtons": [["toImage", "resetScale2d"]],
    "scrollZoom": False,
    "displaylogo": False,
}

COLORS = px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet


class Plotly_plot:
    def __init__(self, df):
        self._df = df
        self._fig = None

    def plotly_plot(self):
        return plotly.offline.plot(
            self._fig, include_plotlyjs=False, output_type="div", config=PLOTLY_CONFIG
        )

    def get_plotly_div(self):
        return self.plotly_plot()


class Insert_plot(Plotly_plot):
    def __init__(self, df):
        super().__init__(df)
        self.set_fig()

    def set_fig(self):
        hist_data = [self._df["size"]]
        group_labels = ["Sample"]
        colors = ["#1F77B4"]

        self._fig = pf.create_distplot(
            hist_data, group_labels, show_hist=False, colors=colors, show_rug=False
        )
        self._fig.update_layout(
            title={
                "text": "Insert size distribution",
                "y": 0.98,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            xaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            xaxis_title="Insert Size(bp)",
            yaxis_title="Density",
            font=dict(size=14, color="Black"),
            showlegend=False,
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF",
            width=470,
            height=313,
        )


class Tss_plot(Plotly_plot):
    def __init__(self, df):
        super().__init__(df)
        self.set_fig()

    def set_fig(self):
        self._fig = px.line(
            self._df,
            x="dists",
            y="agg_tss_scores",
            markers=True,
            color_discrete_sequence=["#1F77B4"],
        )
        self._fig.update_layout(
            title={
                "text": "Aggregate TSS enrichment",
                "y": 0.98,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            xaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            xaxis_title="Distance from TSS(bp)",
            yaxis_title="Aggregate TSS score",
            font=dict(size=14, color="Black"),
            showlegend=False,
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF",
        )


class Peak_plot(Plotly_plot):
    def __init__(self, df):
        super().__init__(df)
        self.set_fig()

    def set_fig(self):
        self._df.sort_values("cell_called", ascending=False, inplace=True)

        self._fig = px.scatter(
            self._df,
            x="fragments",
            y="frac_peak",
            color="cell_called",
            color_discrete_sequence=["grey", "orange"],
        )
        self._fig.update_layout(
            width=470,
            height=313,
            title={
                "text": "Peaks Targeting",
                "y": 0.98,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            xaxis={
                "type": "log",
                "color": "black",
                "gridcolor": "gainsboro",
                "linecolor": "black",
            },
            yaxis={"color": "black", "gridcolor": "gainsboro", "linecolor": "black"},
            xaxis_title="Fragments per Barcode",
            yaxis_title="Fraction Fragments Overlapping Peaks",
            font=dict(size=12, color="Black"),
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF",
            legend=dict(title_text=""),
        )
        self._fig.update_traces(marker_size=3)


class Frag_dis_plot(Plotly_plot):
    def __init__(self, df):
        super().__init__(df)
        self.set_fig()

    @staticmethod
    def get_plot_data(df):
        singlecell_df = df
        fragment_counts = singlecell_df["fragments"].values
        cell_mask = singlecell_df["cell_called"] == "Cells"
        noncell_mask = singlecell_df["cell_called"] == "Non-cells"

        logbinmax = np.ceil(np.log10(fragment_counts.max()))
        xbins = list(
            np.hstack([np.arange(100), np.logspace(np.log10(100), logbinmax, 350)])
        )
        data_subplots = []
        for name, mask, color in zip(
            ["Non-cells", "{} Cells".format("")],
            [noncell_mask, cell_mask],
            ["blue", "orange"],
        ):
            # protect against really low depth samples, mixed up ref, severe under cell calling
            if mask.sum() > 0 and logbinmax > np.log10(100):
                counts, _ = np.histogram(fragment_counts[mask], xbins)
                data_subplots.append(
                    {
                        "name": name,
                        "x": xbins,
                        "y": list(counts),
                        "type": "scatter",
                        "connectgaps": True,
                        "fill": "tozeroy",
                        "line": {"color": color},
                    }
                )

        return data_subplots

    def set_fig(self):
        data_subplots = self.get_plot_data(self._df)

        self._fig = go.Figure()
        self._fig.add_trace(
            go.Scatter(
                x=data_subplots[0]["x"],
                y=data_subplots[0]["y"],
                fill="tozeroy",
                connectgaps=True,
                name="Non-cells",
            )
        )  # fill down to xaxis
        self._fig.add_trace(
            go.Scatter(
                x=data_subplots[1]["x"],
                y=data_subplots[1]["y"],
                fill="tozeroy",
                connectgaps=True,
                name="Cells",
            )
        )  # fill to trace0 y

        self._fig.update_layout(
            width=470,
            height=313,
            title={
                "text": "Fragment Distribution",
                "y": 0.98,
                "x": 0.5,
                "xanchor": "center",
                "yanchor": "top",
            },
            xaxis={
                "type": "log",
                "color": "black",
                "gridcolor": "gainsboro",
                "linecolor": "black",
            },
            yaxis={
                "type": "log",
                "color": "black",
                "gridcolor": "gainsboro",
                "linecolor": "black",
            },
            xaxis_title="Fragments per Barcode",
            yaxis_title="Barcode",
            font=dict(size=12, color="Black"),
            margin=dict(l=50, r=0, t=30, b=30),
            plot_bgcolor="#FFFFFF",
            legend=dict(title_text=""),
        )


class Umap_plot(Plotly_plot):
    def __init__(self, df_umap, feature_name, discrete=True):
        super().__init__(df_umap)

        self.feature_name = feature_name
        self.discrete = discrete
        title_feature_name = feature_name[0].upper() + feature_name[1:]
        self.title = f"UMAP plot Colored by {title_feature_name}"

        self._layout = {}
        self._dot_size = 4
        self._df["size"] = self._dot_size
        self._df["barcode_index"] = list(range(1, len(self._df) + 1))
        self._str_coord1 = "UMAP_1"
        self._str_coord2 = "UMAP_2"
        self.axes_config = {
            "showgrid": True,
            "gridcolor": "#F5F5F5",
            "showline": False,
            "ticks": None,
            "zeroline": True,
            "zerolinecolor": "black",
            "zerolinewidth": 0.7,
        }

        self.scatter_config = {
            "data_frame": df_umap,
            "title": self.title,
            "x": self._str_coord1,
            "y": self._str_coord2,
            "size_max": self._dot_size,
            "hover_data": {
                self._str_coord1: False,
                self._str_coord2: False,
                self.feature_name: True,
                "barcode_index": True,
                "size": False,
            },
            "size": "size",
            "opacity": 0.9,
            "color": self.feature_name,
            "color_discrete_sequence": COLORS,
            "color_continuous_scale": px.colors.sequential.Jet,
        }

    def set_color_scale(self, color_scale):
        self.scatter_config["color_continuous_scale"] = color_scale

    @utils.add_log
    def get_plotly_div(self):
        if self.discrete:
            self.discrete_umap_plot()
        else:
            self.continuous_umap_plot()
        self.update_fig()

        return self.plotly_plot()

    @utils.add_log
    def discrete_umap_plot(self):
        sum_df = self._df.groupby([self.feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = defaultdict(int)
        res_list = []
        for cluster in sorted(self._df[self.feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster] = name
            res_list.append(name)

        self._df[self.feature_name] = self._df[self.feature_name].map(res_dict)

        self._fig = px.scatter(
            **self.scatter_config, category_orders={self.feature_name: res_list}
        )

    @utils.add_log
    def continuous_umap_plot(self):
        self._fig = px.scatter(
            **self.scatter_config,
        )

    def update_fig(self):
        self._fig.update_xaxes(title_text=self._str_coord1, **self.axes_config)

        self._fig.update_yaxes(title_text=self._str_coord2, **self.axes_config)

        self._fig.update_layout(
            self._layout,
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor="#FFFFFF",
            hovermode="closest",
        )


class Tsne_plot(Plotly_plot):
    def __init__(self, df_tsne, feature_name, discrete=True):
        super().__init__(df_tsne)

        self.feature_name = feature_name
        self.discrete = discrete
        title_feature_name = feature_name[0].upper() + feature_name[1:]
        self.title = f"t-SNE plot Colored by {title_feature_name}"

        self._layout = {}
        self._dot_size = 4
        self._df["size"] = self._dot_size
        self._df["barcode_index"] = list(range(1, len(self._df) + 1))
        self._str_coord1 = "tSNE_1"
        self._str_coord2 = "tSNE_2"
        self.axes_config = {
            "showgrid": True,
            "gridcolor": "#F5F5F5",
            "showline": False,
            "ticks": None,
            "zeroline": True,
            "zerolinecolor": "black",
            "zerolinewidth": 0.7,
        }

        self.scatter_config = {
            "data_frame": df_tsne,
            "title": self.title,
            "x": self._str_coord1,
            "y": self._str_coord2,
            "size_max": self._dot_size,
            "hover_data": {
                self._str_coord1: False,
                self._str_coord2: False,
                self.feature_name: True,
                "barcode_index": True,
                "size": False,
            },
            "size": "size",
            "opacity": 0.9,
            "color": self.feature_name,
            "color_discrete_sequence": COLORS,
            "color_continuous_scale": px.colors.sequential.Jet,
        }

    def set_color_scale(self, color_scale):
        self.scatter_config["color_continuous_scale"] = color_scale

    @utils.add_log
    def get_plotly_div(self):
        if self.discrete:
            self.discrete_tsne_plot()
        else:
            self.continuous_tsne_plot()
        self.update_fig()

        return self.plotly_plot()

    @utils.add_log
    def discrete_tsne_plot(self):
        sum_df = self._df.groupby([self.feature_name]).agg("count").iloc[:, 0]
        percent_df = sum_df.transform(lambda x: round(x / sum(x) * 100, 2))
        res_dict = dict()
        res_list = []
        for cluster in sorted(self._df[self.feature_name].unique()):
            name = f"{cluster}({percent_df[cluster]}%)"
            res_dict[cluster] = name
            res_list.append(name)

        self._df[self.feature_name] = self._df[self.feature_name].map(res_dict)

        self._fig = px.scatter(
            **self.scatter_config, category_orders={self.feature_name: res_list}
        )

    @utils.add_log
    def continuous_tsne_plot(self):
        self._fig = px.scatter(
            **self.scatter_config,
        )

    def update_fig(self):
        self._fig.update_xaxes(title_text=self._str_coord1, **self.axes_config)

        self._fig.update_yaxes(title_text=self._str_coord2, **self.axes_config)

        self._fig.update_layout(
            self._layout,
            title={"text": self.title, "x": 0.5, "y": 0.95, "font": {"size": 15}},
            plot_bgcolor="#FFFFFF",
            hovermode="closest",
        )
