import pathlib
import pandas as pd
from plotnine import (
    ggplot,
    aes,
    geom_tile,
    labs,
    theme,
    element_text,
    geom_text,
    theme_bw,
    scale_color_manual,
    facet_grid,
    element_blank,
)
from k4neo.annotator import TUMOR_TISSUE, IMMUNO_PRIVILIGED_TISSUE


class Plotter:
    """Generates a plot from a tidy sample rate DataFrame generated with k4neo"""

    def __init__(
        self,
        healthy_sample_rate: str,
        tumor_sample_rate: str,
        row: str = "junc_id",
        column: str = "tissue",
        value: str = "sample_rate",
    ):
        """
        Initializes Plotter class.

        Args:
            df (pd.DataFrame): Input DataFrame in tidy format.
            row (str): Column to use for the y-axis (rows).
            column (str): Column to use for the x-axis (columns).
            value (str): Column to use for the tile fill values.
        """
        self.healthy_df = pd.read_csv(healthy_sample_rate, sep="\t")
        self.healthy_df = self.healthy_df[self.healthy_df["developmental_stage"] == "adult"]
        self.tumor_df = pd.read_csv(tumor_sample_rate, sep="\t")
        self.row = row
        self.column = column
        self.value = value

    def generate_plot_df(self):
        """
        Get necessary columns from dataframe and merge into one df for printing
        """
        df_healthy = self.healthy_df[[self.row, self.column, self.value]]
        df_healthy["tissue_group"] = df_healthy.apply(
            lambda row: (
                "immuno-\nprivileged" if row.tissue in IMMUNO_PRIVILIGED_TISSUE else "healthy"
            ),
            axis=1,
        )

        df_cancer = self.tumor_df[[self.row, self.column, self.value, "disease"]]
        df_cancer = df_cancer[df_cancer["disease"].isin(TUMOR_TISSUE)]
        df_cancer.rename(columns={"disease": "tissue_group"}, inplace=True)

        return pd.concat([df_healthy, df_cancer])

    def plot(self, output_path: str = None, dpi: int = 300, width=12, height=9):
        """
        Generates and optionally saves a heatmap plot.

        Args:
            output_path (str): Path to save the plot (e.g., PNG, PDF). If None, plot is shown.
            dpi (int): Resolution of the saved image.
            title (str): Title of the plot.

        Returns:
            ggplot: The generated plotnine ggplot object.
        """
        df = self.generate_plot_df()
        df["p_group"] = pd.cut(
            df["sample_rate"], (0, 0.35, 1), labels=("low", "high"), include_lowest=True
        )
        plot = (
            ggplot(df, aes(self.column, self.row, fill=self.value))
            + geom_tile(aes(width=0.95, height=0.95))
            + geom_text(
                aes(label=f"round(sample_rate,2)", color="p_group"), size=5, show_legend=False
            )
            + theme_bw()
            + scale_color_manual(["white", "black"])
            + facet_grid(
                cols="tissue_group",
                scales="free",
                space="free",
            )
            + theme(
                axis_text_x=element_text(rotation=45, hjust=1),
                axis_title_y=element_blank(),
                axis_ticks_y=element_blank(),
                panel_grid_major=element_blank(),
                panel_grid_minor=element_blank(),
                legend_position="right",
            )
        )

        if output_path:
            plot.save(output_path, dpi=dpi, width=width, height=height)

        return plot
