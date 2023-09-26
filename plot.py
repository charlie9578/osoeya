import matplotlib as mpl

from pyproj import Transformer
from bokeh.models import WMTSTileSource, ColumnDataSource
from bokeh.palettes import Category10, viridis
from bokeh.plotting import figure


def luminance(rgb):
    """Calculates the brightness of an rgb 255 color. See https://en.wikipedia.org/wiki/Relative_luminance

    Args:
        rgb(:obj:`tuple`): 255 (red, green, blue) tuple

    Returns:
        luminance(:obj:`scalar`): relative luminance

    Example:

        .. code-block:: python

            >>> rgb = (255,127,0)
            >>> luminance(rgb)
            0.5687976470588235

            >>> luminance((0,50,255))
            0.21243529411764706

    """

    luminance = (0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]) / 255

    return luminance


def color_to_rgb(color):
    """Converts named colors, hex and normalised RGB to 255 RGB values

    Args:
        color(:obj:`color`): RGB, HEX or named color

    Returns:
        rgb(:obj:`tuple`): 255 RGB values

    Example:

        .. code-block:: python

            >>> color_to_rgb("Red")
            (255, 0, 0)

            >>> color_to_rgb((1,1,0))
            (255,255,0)

            >>> color_to_rgb("#ff00ff")
            (255,0,255)
    """

    if isinstance(color, tuple):
        if max(color) > 1:
            color = tuple([i / 255 for i in color])

    rgb = mpl.colors.to_rgb(color)

    rgb = tuple([int(i * 255) for i in rgb])

    return rgb


def plot_assets(
    asset_df,
    selected_asset=None,
    tile_name="OpenMap",
    plot_width=800,
    plot_height=800,
    marker_size=14,
    kwargs_for_figure={},
    kwargs_for_marker={},
):
    """Plot the windfarm spatially on a map using the Bokeh plotting libaray.

    Args:
        asset_df(:obj:`pd.DataFrame`): PlantData.asset object containing the asset metadata.
        tile_name(:obj:`str`): tile set to be used for the underlay, e.g. OpenMap, ESRI, OpenTopoMap
        plot_width(:obj:`scalar`): width of plot
        plot_height(:obj:`scalar`): height of plot
        marker_size(:obj:`scalar`): size of markers
        kwargs_for_figure(:obj:`dict`): additional figure options for advanced users, see Bokeh docs
        kwargs_for_marker(:obj:`dict`): additional marker options for advanced users, see Bokeh docs. We have some custom behavior around the "fill_color" attribute. If "fill_color" is not defined, OpenOA will use an internally defined color pallete. If "fill_color" is the name of a column in the asset table, OpenOA will use the value of that column as the marker color. Otherwise, "fill_color" is passed through to Bokeh.

    Returns:
        Bokeh_plot(:obj:`axes handle`): windfarm map

    """
    """
    TODO: Add this back in when it can be debugged
    lats, lons = transformer.transform(...) -> ValueError: not enough values to unpack (expected 2, got 0)

    Example:
        .. bokeh-plot::

            import pandas as pd
            from bokeh.plotting import figure, output_file, show

            from openoa.utils.plot import plot_windfarm

            from examples import project_ENGIE

            # Load plant object
            project = project_ENGIE.prepare("../examples/data/la_haute_borne")

            # Create the bokeh wind farm plot
            show(plot_windfarm(project.asset, tile_name="ESRI", plot_width=600, plot_height=600))
    """

    # See https://wiki.openstreetmap.org/wiki/Tile_servers for various tile services
    MAP_TILES = {
        "OpenMap": WMTSTileSource(url="http://c.tile.openstreetmap.org/{Z}/{X}/{Y}.png"),
        "ESRI": WMTSTileSource(
            url="https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{Z}/{Y}/{X}.jpg"
        ),
        "OpenTopoMap": WMTSTileSource(url="https://tile.opentopomap.org/{Z}/{X}/{Y}.png"),
    }

    # Use pyproj to transform longitude and latitude into web-mercator and add to a copy of the asset dataframe
    TRANSFORM_4326_TO_3857 = Transformer.from_crs("EPSG:4326", "EPSG:3857")

    asset_df["x"], asset_df["y"] = TRANSFORM_4326_TO_3857.transform(
        asset_df["latitude"], asset_df["longitude"]
    )
    asset_df["coordinates"] = tuple(zip(asset_df["latitude"], asset_df["longitude"]))

    # Define default and then update figure and marker options based on kwargs
    figure_options = {
        "tools": "save,hover,pan,wheel_zoom,reset,help",
        "x_axis_label": "Longitude",
        "y_axis_label": "Latitude",
        "match_aspect": True,
        "tooltips": [("plantName", "@plantName"), 
                     ("technology", "@technology"), 
                     ("capacity MW", "@nameplateCapacityMW"), 
                     ("COD", "@operatingYearMonth"),
                     ("status", "@statusDescription"),
                     ("(Lat,Lon)", "@coordinates")],
    }
    figure_options.update(kwargs_for_figure)

    marker_options = {
        "marker": "circle_y",
        "line_width": 1,
        "alpha": 0.8,
        "fill_color": "auto_fill_color",
        "line_color": "auto_line_color",
        "legend_group": "technology",
    }
    marker_options.update(kwargs_for_marker)

    # If an asset is selected, re-assign its legend
    if selected_asset is not None:
        asset_df.loc[asset_df["plantid"]==selected_asset["plantid"],"technology"] = "Selected plant"

    # Create an appropriate fill color map and contrasting line color
    if marker_options["fill_color"] == "auto_fill_color":
        color_grouping = marker_options["legend_group"]

        asset_df = asset_df.sort_values(color_grouping)

        if len(set(asset_df[color_grouping])) <= 10:
            color_palette = list(Category10[10])
        else:
            color_palette = viridis(len(set(asset_df[color_grouping])))

        color_mapping = dict(zip(set(asset_df[color_grouping]), color_palette))
        asset_df["auto_fill_color"] = asset_df[color_grouping].map(color_mapping)
        asset_df["auto_fill_color"] = asset_df["auto_fill_color"].apply(color_to_rgb)
        asset_df["auto_line_color"] = [
            "black" if luminance(color) > 0.5 else "white" for color in asset_df["auto_fill_color"]
        ]

    else:
        if marker_options["fill_color"] in asset_df.columns:
            asset_df[marker_options["fill_color"]] = asset_df[marker_options["fill_color"]].apply(
                color_to_rgb
            )
            asset_df["auto_line_color"] = [
                "black" if luminance(color) > 0.5 else "white"
                for color in asset_df[marker_options["fill_color"]]
            ]

        else:
            asset_df["auto_line_color"] = "black"

    # Create the bokeh data source
    source = ColumnDataSource(asset_df)

    # Create a bokeh figure with tiles
    plot_map = figure(
        width=plot_width,
        height=plot_height,
        x_axis_type="mercator",
        y_axis_type="mercator",
        **figure_options,
    )

    plot_map.add_tile(MAP_TILES[tile_name])

    # Plot the asset devices
    plot_map.scatter(x="x", y="y", source=source, size=marker_size, **marker_options)

    # Zoom to selected asset
    if selected_asset is not None:
        x_zoom, y_zoom = TRANSFORM_4326_TO_3857.transform(
            [selected_asset["latitude"]], [selected_asset["longitude"]]
        )

        plot_map.x_range.start = x_zoom[0] - 100000
        plot_map.x_range.end = x_zoom[0] + 100000
        plot_map.y_range.start = y_zoom[0] - 100000
        plot_map.y_range.end = y_zoom[0] + 100000
        

    return plot_map
