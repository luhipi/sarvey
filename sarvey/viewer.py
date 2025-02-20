#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2025 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
#
# This software was developed together with FERN.Lab (fernlab@gfz-potsdam.de) in the context
# of the SAR4Infra project with funds of the German Federal Ministry for Digital and
# Transport and contributions from Landesamt fuer Vermessung und Geoinformation
# Schleswig-Holstein and Landesbetrieb Strassenbau und Verkehr Schleswig-Holstein.
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# Important: This package uses PyMaxFlow. The core of PyMaxflows library is the C++
# implementation by Vladimir Kolmogorov. It is also licensed under the GPL, but it REQUIRES that you
# cite [BOYKOV04] (see LICENSE) in any resulting publication if you use this code for research purposes.
# This requirement extends to SARvey.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

"""Viewer Module for SARvey."""
import os
from typing import Any, Optional
from logging import Logger
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib.collections import PathCollection
from matplotlib import widgets
from matplotlib.backend_bases import MouseButton
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import numpy as np
from scipy.spatial import KDTree
import datetime
import cmcrameri as cmc

from mintpy.objects.colors import ColormapExt  # for DEM_print colormap
from mintpy.utils import readfile
from mintpy.utils.plot import auto_flip_direction
from miaplpy.objects.slcStack import slcStack

from sarvey.objects import AmplitudeImage, Points, BaseStack
import sarvey.utils as ut





def plotIfgs(*, phase: np.ndarray, coord: np.ndarray, spatial_ref_idx: int = None, ttl: str = None,
             cmap: str = "romaO"):
    """Plot one interferogram per subplot.

    Parameters
    ----------
    phase: np.ndarray
        phase per point and ifg, e.g. wrapped or unwrapped phase (dim: no. psPoints x no. ifgs)
    coord: np.ndarray
        coordinates of the psPoints, e.g. pixel or lat lon (dim: no. psPoints x 2)
    spatial_ref_idx: int
        index of the spatial reference point (default: None)
    ttl: str
        title for the figure (default: None)
    cmap: str
        colormap name (default: "romaO")
    """
    num_ifgs = phase.shape[1]
    min_val = np.min(phase)
    max_val = np.max(phase)
    fig, axs = plt.subplots(np.ceil(np.sqrt(num_ifgs + 1)).astype(np.int32),
                            np.ceil(np.sqrt(num_ifgs + 1)).astype(np.int32))
    sc = None
    for i, ax in enumerate(axs.flat):
        if i < num_ifgs:
            sc = ax.scatter(coord[:, 1], coord[:, 0], c=phase[:, i],
                            vmin=min_val, vmax=max_val, s=1, cmap=cmc.cm.cmaps[cmap])
            ax.axes.set_xticks([])
            ax.axes.set_yticks([])
            if spatial_ref_idx is not None:
                ax.plot(coord[spatial_ref_idx, 1],
                        coord[spatial_ref_idx, 0], 'k*')
        elif i == num_ifgs:
            plt.colorbar(sc, cax=ax)
        else:
            ax.set_visible(False)
    if ttl is not None:
        fig.suptitle(ttl)


def plotScatter(*, value: np.ndarray, coord: np.ndarray, bmap_obj: AmplitudeImage = None, ttl: str = None,
                unit: str = None, s: float = 5.0, cmap: str = "batlow", symmetric: bool = False,
                logger: Logger, **kwargs: Any):
    """Plot a scatter map for given value.

    Parameters
    ----------
    value: np.ndarray
        value to be plotted per point giving the colour of the point (dim: no. points x 1)
    coord: np.ndarray
        coordinates of the points, e.g. radar or lat lon (dim: no. points x 2). If bmapObj is given,
        the coordinates must be radar coordinates!
    bmap_obj: AmplitudeImage
        instance of amplitudeImage for plotting background image (default: None)
    ttl: str
        title for the figure (default: None)
    unit: str
        unit as title for the colorbar axis (default: None)
    s: float
        size of the scatter points (default: 5.0)
    cmap: str
        colormap (default: "batlow")
    symmetric: bool
        plot symmetric colormap extend, i.e. abs(vmin) == abs(vmax) (default: False)
    logger: Logger
        logging Handler
    kwargs: Any
        additional keyword arguments for scatter plot

    Returns
    -------
    fig: plt.Figure
        current figure,
    ax: plt.Axes
        current axis
    cb: plt.colorbar
        current colorbar
    """
    if bmap_obj is not None:
        ax = bmap_obj.plot(logger=logger)
        fig = plt.gcf()
    else:
        fig = plt.figure()
        ax = fig.add_subplot()

    if symmetric:
        v_range = np.max(np.abs(value.ravel()))
        sc = ax.scatter(coord[:, 1], coord[:, 0], c=value, s=s, cmap=cmc.cm.cmaps[cmap],
                        vmin=-v_range, vmax=v_range)
    else:
        sc = ax.scatter(coord[:, 1], coord[:, 0], c=value, s=s, cmap=cmc.cm.cmaps[cmap], **kwargs)

    cb = plt.colorbar(sc, ax=ax, pad=0.03, shrink=0.5)
    cb.ax.set_title(unit)
    ax.set_title(ttl)
    plt.tight_layout()
    return fig, ax, cb


def plotColoredPointNetwork(*, x: np.ndarray, y: np.ndarray, arcs: np.ndarray, val: np.ndarray, ax: plt.Axes = None,
                            linewidth: float = 2, cmap: str = "vik", clim: tuple = None):
    """Plot a network of points with colored arcs.

    Parameters
    ----------
    x: np.ndarray
        x-coordinates of the points (dim: no. points x 1)
    y: np.ndarray
        y-coordinates of the points (dim: no. points x 1)
    arcs: np.ndarray
        indices of the points to be connected (dim: no. arcs x 2)
    val: np.ndarray
        values for the color of the arcs (dim: no. arcs x 1)
    ax: plt.Axes
        axis for plotting (default: None)
    linewidth: float
        line width of the arcs (default: 2)
    cmap: str
        name of the colormap (default: "vik")
    clim: tuple
        color limits for the colormap (default: None)

    Returns
    -------
    ax: plt.Axes
        current axis
    cbar: plt.colorbar
        current colorbar
    """
    if ax is None:
        fig = plt.figure(figsize=(15, 5))
        ax = fig.add_subplot()
    else:
        fig = ax.get_figure()
    ax.scatter(x, y, s=3.5, c=np.ones_like(x))

    if clim is None:
        norm = Normalize(vmin=min(val), vmax=max(val))
    else:
        norm = Normalize(vmin=clim[0], vmax=clim[1])

    mapper = cm.ScalarMappable(norm=norm, cmap=cmc.cm.cmaps[cmap])
    mapper_list = [mapper.to_rgba(v) for v in val]
    for m in range(arcs.shape[0]):
        x_val = [x[arcs[m, 0]], x[arcs[m, 1]]]
        y_val = [y[arcs[m, 0]], y[arcs[m, 1]]]

        ax.plot(x_val, y_val, linewidth=linewidth, c=mapper_list[m])
    cbar = fig.colorbar(mapper, ax=ax, pad=0.03, shrink=0.5)

    return ax, cbar


def plotGridFromBoxList(*, box_list: list, ax: plt.Axes = None, edgecolor: str = "k", linewidth: float = 1):
    """Plot a grid into an axis.

    Parameters
    ----------
    box_list: list
        boxes to be plotted. box_list can be created with 'splitImageIntoBoxesRngAz' or 'splitImageIntoBoxes'
    ax: plt.Axes
        axis for plotting (default: None)
    edgecolor: str
        edge color for the boxes (default: "k")
    linewidth: float
        line width for the boxes (default: 1)

    Returns
    -------
    ax: plt.Axes
        current axis
    """
    if ax is None:
        fig = plt.figure()
        ax = fig.add_subplot()

    for box in box_list:
        rect = patches.Rectangle((box[0], box[1]), box[2] - box[0], box[3] - box[1], linewidth=linewidth,
                                 edgecolor=edgecolor, facecolor="none")
        ax.add_patch(rect)
    return ax


class LineSelector:
    """LineSelector."""

    def __init__(self, ax):
        """Init.

        Parameters
        ----------
        ax: plt.Axes
            axis
        """
        self.ax = ax
        self.unselected_color = 'k'
        self.points = []
        self.first_point = None
        self.second_point = None
        self.select_amplitude = True
        self.line = None
        self.line_point_indices = [None, None]

    def plotImageAcquisitions(self, slc_stack_obj: slcStack, date_list: list):
        """Initialize figure by plotting the baseline plot of the image acquisitions.

        Parameters
        ----------
        slc_stack_obj: slcStack
            instance of slcStack
        date_list: list
            list of dates
        """
        for date_idx in range(slc_stack_obj.numDate):
            point = self.ax.plot(date_list[date_idx], slc_stack_obj.pbase[date_idx], "ko")[0]
            self.points.append(point)
        self.ax.set_title("Image acquistions")
        self.ax.set_xlabel("Date")
        self.ax.set_ylabel("Perpendicular baseline [m]")

    def onClick(self, event):
        """Select amplitude or interferogram.

        Parameters
        ----------
        event: MouseEvent
            mouse event
        """
        if event.inaxes != self.ax:
            return None, None

        if self.select_amplitude:
            selected_date_idx = None
            for date_idx in range(len(self.points)):
                point = self.points[date_idx]
                if point.contains(event)[0]:
                    point.set_color('r')
                    selected_date_idx = date_idx
                else:
                    point.set_color(self.unselected_color)
            plt.draw()
            return selected_date_idx, None

        else:  # select interferograms
            for date_idx in range(len(self.points)):
                point = self.points[date_idx]
                if point.contains(event)[0]:
                    if event.button == 1:  # left click
                        self.first_point = point
                        self.line_point_indices[0] = date_idx
                    if event.button == 3:  # right click
                        self.second_point = point
                        self.line_point_indices[1] = date_idx
                    break

            if (self.line_point_indices[0] is not None) and (self.line_point_indices[1] is not None):
                line_coord = [[self.first_point.get_xdata(), self.second_point.get_xdata()],
                              [self.first_point.get_ydata(), self.second_point.get_ydata()]]

                if self.line is None:  # only for first plot
                    self.line = Line2D(line_coord[0], line_coord[1], color='b')
                    self.ax.add_line(self.line)
                else:
                    self.line.set_data(line_coord[0], line_coord[1])
                plt.draw()
            else:
                print("First image is selected. Now select the second image.")

            return self.line_point_indices[0], self.line_point_indices[1]

    def onCheck(self, label):
        """Check if amplitude or interferogram is selected.

        Parameters
        ----------
        label: str
            label of the radio button
        """
        if label == 'Amplitude':
            self.select_amplitude = True
            # remove line
            if self.line is not None:
                self.line.remove()
                self.line = None
                self.line_point_indices = [None, None]
                plt.draw()
        else:
            self.select_amplitude = False
            print("Create interferograms by clicking on acquisitions."
                  "\nFirst acquistion: LEFT mouse click"
                  "\nSecond acquistion: RIGHT mouse click\n")
            print("Start with selecting the first image.")
            for point in self.points:
                point.set_color(self.unselected_color)


class ImageViewer:
    """ImageViewer."""

    def __init__(self, slc_stack_obj: slcStack, line_selector: LineSelector):
        """Init.

        Parameters
        ----------
        slc_stack_obj: slcStack
            instance of slcStack
        line_selector: LineSelector
            instance of LineSelector
        """
        self.slc_stack_obj = slc_stack_obj
        self.line_selector = line_selector

    def initializeImage(self):
        """InitializeImage."""
        self.fig_img = plt.figure()
        self.ax_img = plt.subplot()
        self.img = self.ax_img.imshow(np.zeros((self.slc_stack_obj.length, self.slc_stack_obj.width), dtype=np.float32),
                                      interpolation="nearest")
        auto_flip_direction(self.slc_stack_obj.metadata, ax=self.ax_img, print_msg=True)
        self.cb = plt.colorbar(self.img, ax=self.ax_img, pad=0.03, shrink=0.5)
        self.ax_img.set_xlabel("Range")
        self.ax_img.set_ylabel("Azimuth")

    def plotImage(self, event):
        """Plot either amplitude image or interferogram.

        Parameters
        ----------
        event: MouseEvent
            mouse event
        """
        date1_idx, date2_idx = self.line_selector.onClick(event)
        if date1_idx is None:
            return
        if self.line_selector.select_amplitude:
            slc = self.slc_stack_obj.read(datasetName=self.slc_stack_obj.dateList[date1_idx], print_msg=False)
            ampl = np.abs(slc)
            ampl[ampl == 0] = np.nan
            ampl = np.log10(ampl)
            self.img.set_data(ampl)
            self.img.set_clim(np.nanmin(ampl), np.nanmax(ampl))
            self.img.set_cmap('gray')
            d = self.slc_stack_obj.dateList[date1_idx]
            d = datetime.date(year=int(d[:4]), month=int(d[4:6]), day=int(d[6:]))
            self.ax_img.set_title(f"{d}")
        else:
            if (date1_idx is None) or (date2_idx is None):
                return
            slc1 = self.slc_stack_obj.read(datasetName=self.slc_stack_obj.dateList[date1_idx], print_msg=False)
            slc2 = self.slc_stack_obj.read(datasetName=self.slc_stack_obj.dateList[date2_idx], print_msg=False)
            ifg = slc2 * np.conjugate(slc1)
            self.img.set_data(np.angle(ifg))
            self.img.set_clim(-np.pi, np.pi)
            self.img.set_cmap(cmc.cm.cmaps["romaO"])

            d1 = self.slc_stack_obj.dateList[date1_idx]
            d2 = self.slc_stack_obj.dateList[date2_idx]
            d1 = datetime.date(year=int(d1[:4]), month=int(d1[4:6]), day=int(d1[6:]))
            d2 = datetime.date(year=int(d2[:4]), month=int(d2[4:6]), day=int(d2[6:]))
            days = self.slc_stack_obj.tbase[date2_idx] - self.slc_stack_obj.tbase[date1_idx]
            perp_base = np.abs(self.slc_stack_obj.pbase[date2_idx] - self.slc_stack_obj.pbase[date1_idx])
            self.ax_img.set_title(f"Interferogram:\n{d1}   -   {d2}\n{int(days)} days, {int(perp_base)} m")
        self.fig_img.canvas.draw()


class TimeSeriesViewer:
    """TimeSeriesViewer."""

    def __init__(self, *, point_obj: Points, vel_scale: str = "mm", input_path: str, logger: Logger):
        """Init."""
        self.sc = None
        self.point_obj = point_obj
        self.ts_point_marker = None  # for ts point marker
        self.ts_point_idx = 0  # index of ts_point
        self.ts_refpoint_marker = None  # for reference point marker
        self.logger = logger
        self.ts_refpoint_idx = None  # index of reference point
        self.vel_scale = vel_scale
        scale_dict = {"mm": 1000, "cm": 100, "dm": 10, "m": 1}
        if self.vel_scale not in scale_dict.keys():
            raise ValueError(f"Invalid argument: '{self.vel_scale}'")
        self.scale = scale_dict[self.vel_scale]
        self.tree = KDTree(self.point_obj.coord_xy)
        self.tree_coord_map = KDTree(self.point_obj.coord_map)
        self.mask_ccs = np.isnan(self.point_obj.phase).sum(axis=1) == 0
        self.tree_ccs = KDTree(self.point_obj.coord_xy[self.mask_ccs])
        if point_obj.ifg_net_obj.dates is not None:
            self.times = [datetime.date.fromisoformat(date) for date in point_obj.ifg_net_obj.dates]
        else:  # backwards compatible, if ifg_net_obj does not contain dates
            self.times = point_obj.ifg_net_obj.tbase

        vel, demerr, ref_atmo, coherence, omega, v_hat = ut.estimateParameters(obj=self.point_obj, ifg_space=False)
        self.vel = vel
        self.demerr = demerr
        self.ref_atmo = ref_atmo

        self.bmap_obj = AmplitudeImage(file_path=os.path.join(os.path.dirname(self.point_obj.file_path),
                                                              "background_map.h5"))
        self.bmap_obj.open()
        self.height = readfile.read(os.path.join(input_path, "geometryRadar.h5"), datasetName='height')[0]

        temp_coh_obj = BaseStack(
            file=os.path.join(os.path.dirname(self.point_obj.file_path), "temporal_coherence.h5"),
            logger=logger)
        self.temp_coh_img = temp_coh_obj.read(dataset_name="temp_coh")

        self.font_size = 10
        plt.rc('font', size=self.font_size)  # controls default text size
        plt.rc('axes', titlesize=self.font_size)  # fontsize of the title
        plt.rc('axes', labelsize=self.font_size)  # fontsize of the x and y labels
        plt.rc('xtick', labelsize=self.font_size)  # fontsize of the x tick labels
        plt.rc('ytick', labelsize=self.font_size)  # fontsize of the y tick labels
        plt.rc('legend', fontsize=self.font_size)  # fontsize of the legend

        self.initFigureMap()
        self.initFigureTimeseries()
        self.plotMap(val=None)
        self.fig1.canvas.mpl_connect('button_press_event', self.onClick)
        plt.show()

    def initFigureMap(self):
        """InitFigureMap."""
        self.fig1 = plt.figure()
        self.ax_img = self.fig1.subplots(1, 1)

        self.ax_cb = self.fig1.add_axes((0.93, 0.6, 0.015, 0.15))  # (left, bottom, width, height)
        self.cb = self.fig1.colorbar(self.sc,
                                     cax=self.ax_cb,
                                     ax=self.ax_img,
                                     pad=0.03,
                                     shrink=0.8,
                                     aspect=10,
                                     orientation='vertical')

        # slider for TCS time span
        self.ax_date_slider = self.fig1.add_axes((0.125, 0.05, 0.75, 0.03))  # (left, bottom, width, height)
        self.date_slider = widgets.RangeSlider(
            ax=self.ax_date_slider,
            label='Image Range',
            valmin=0,
            valmax=self.point_obj.ifg_net_obj.num_images - 1,
            valinit=(0, self.point_obj.ifg_net_obj.num_images - 1),
            valstep=1,
        )

        def label_formatter(val):
            """Map slider value to string."""
            return f"{self.times[int(val[0])]}\n{self.times[int(val[1])]}"

        self.date_slider.valtext.set_text(label_formatter(self.date_slider.val))  # Set initial text
        self.date_slider.on_changed(lambda val: self.date_slider.valtext.set_text(label_formatter(val)))
        self.date_slider.on_changed(self.plotMap)
        self.date_slider.on_changed(self.updateBackgroundTimeseriesPlot)

        # add button to select reference point
        self.set_reference_point = False
        self.ax_button = self.fig1.add_axes((0.05, 0.9, 0.1, 0.08))  # (left, bottom, width, height)
        self.button_mask = widgets.Button(ax=self.ax_button, label='Select\nReference', image=None, color='1')
        self.button_mask.on_clicked(self.updateButtonStatus)

        # add radiobutton to select parameter
        self.ax_radio_par = self.fig1.add_axes((0.15, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.rb_par = widgets.RadioButtons(self.ax_radio_par, labels=['Velocity', 'DEM correction', 'None'], active=0)
        self.rb_par.on_clicked(self.plotMap)

        # add radiobutton to select background image
        self.ax_radio_backgr = self.fig1.add_axes((0.35, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.rb_backgr = widgets.RadioButtons(self.ax_radio_backgr, labels=['Amplitude', 'DEM', 'Coherence', 'None'],
                                              active=0)
        self.rb_backgr.on_clicked(self.plotMap)

        # add info box with info about velocity and DEM error of selected pixel
        self.ax_info_box = self.fig1.add_axes((0.55, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.text_obj_time = self.ax_info_box.text(0.1, 0.1, "")
        self.ax_info_box.set_xticks([], [])
        self.ax_info_box.set_yticks([], [])

        # add variable for axis of slider controlling the visualized coherence background image
        self.ax_slide_coh = None
        self.sl_last_val = 0.0
        self.sl_coh = None

        # add neighbourhood markers
        self.neighb_markers = None
        self.neighb_ts_lines = list()
        self.ax_radius = self.fig1.add_axes((0.75, 0.9, 0.1, 0.08))  # (left, bottom, width, height)
        self.txt_radius = widgets.TextBox(self.ax_radius, 'Radius [m]', label_pad=-1, initial="30")
        self.txt_radius.on_submit(self.updateNeighbourhood)
        self.txt_radius.on_submit(self.plotPointTimeseries)

        # add radiobutton to show location of CCS
        self.ax_radio_ccs = self.fig1.add_axes((0.85, 0.9, 0.1, 0.08))  # (left, bottom, width, height)
        self.rb_ccs = widgets.RadioButtons(ax=self.ax_radio_ccs, labels=['None', 'CCS'], active=0)
        self.rb_ccs.on_clicked(self.plotMap)

    def initFigureTimeseries(self):
        """InitFigureTimeseries."""
        self.fig2 = plt.figure(figsize=(15, 5))
        self.ax_ts = self.fig2.subplots(1, 1)

        # add radiobutton for fitting linear model
        self.ax_radio_fit = self.fig2.add_axes((0.125, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.rb_fit = widgets.RadioButtons(self.ax_radio_fit, labels=['None', 'Linear fit'], active=0)

        # add radiobutton for selecting baseline type
        self.ax_radio_baselines = self.fig2.add_axes((0.325, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.rb_baselines = widgets.RadioButtons(
            self.ax_radio_baselines,
            labels=['Temporal baseline', 'Perpendicular baseline'],
            active=0
        )

        # add check box for removing phase due to parameters
        self.ax_cbox_par = self.fig2.add_axes((0.525, 0.9, 0.2, 0.08))  # (left, bottom, width, height)
        self.cbox_par = widgets.CheckButtons(
            self.ax_cbox_par,
            ["Velocity", "DEM correction"],
            actives=[True, False]
        )
        self.rb_fit.on_clicked(self.plotPointTimeseries)
        self.rb_baselines.on_clicked(self.plotPointTimeseries)
        self.cbox_par.on_clicked(self.plotPointTimeseries)

        # initialize time series
        self.ts_line = self.ax_ts.plot(self.times, np.zeros_like(self.times), '.')[0]
        self.lin_fit_line = self.ax_ts.plot([], [], '-k')[0]

        # add highlighted area which indicates the selected time span
        self.highlight_time_span = self.ax_ts.axvspan(
            xmin=self.times[0],
            xmax=self.times[-1],
            color='grey',
            alpha=0.1
        )
        self.ax_ts.set_xlim(self.times[0], self.times[-1])

    def plotMap(self, val: object):  # val seems to be unused, but its necessary for the function to work.
        """Plot velocity map and time series."""
        flag_initial_plot = (0.0, 1.0) == self.ax_img.get_xlim()
        ax_img_xlim = None
        ax_img_ylim = None
        if not flag_initial_plot:
            ax_img_xlim = self.ax_img.get_xlim()
            ax_img_ylim = self.ax_img.get_ylim()

        self.ax_img.cla()

        # get selected background from radiobutton
        if self.rb_backgr.value_selected == "Amplitude":
            self.ax_img = self.bmap_obj.plot(ax=self.ax_img, logger=self.logger)
            if self.ax_slide_coh is not None:
                self.sl_last_val = self.sl_coh.val
                self.ax_slide_coh.remove()
                self.ax_slide_coh = None
        if self.rb_backgr.value_selected == "DEM":
            self.ax_img.imshow(self.height, cmap=ColormapExt('DEM_print').colormap)
            meta = {"ORBIT_DIRECTION": self.bmap_obj.orbit_direction}
            auto_flip_direction(meta, ax=self.ax_img, print_msg=False)
            self.ax_img.set_xlabel("Range")
            self.ax_img.set_ylabel("Azimuth")
            if self.ax_slide_coh is not None:
                self.sl_last_val = self.sl_coh.val
                self.ax_slide_coh.remove()
                self.ax_slide_coh = None
        if self.rb_backgr.value_selected == "Coherence":
            if self.ax_slide_coh is None:
                # add slider to change value of coherence for background map
                self.ax_slide_coh = self.fig1.add_axes((0.425, 0.85, 0.2, 0.03))  # (left, bottom, width, height)
                self.sl_coh = widgets.Slider(self.ax_slide_coh,
                                             label='Coherence',
                                             valmin=0.0,
                                             valmax=1.0,
                                             valinit=self.sl_last_val,
                                             valfmt="%.1f")

            self.ax_img.imshow(self.temp_coh_img,
                               cmap=cmc.cm.cmaps["grayC"],
                               vmin=np.round(self.sl_coh.val, decimals=1),
                               vmax=1)
            meta = {"ORBIT_DIRECTION": self.bmap_obj.orbit_direction}
            auto_flip_direction(meta, ax=self.ax_img, print_msg=False)
            self.ax_img.set_xlabel("Range")
            self.ax_img.set_ylabel("Azimuth")
        if self.rb_backgr.value_selected == "None":
            self.ax_img.imshow(np.ones_like(self.height, dtype=np.int8), cmap=cmc.cm.cmaps["grayC"], vmin=0, vmax=1)
            meta = {"ORBIT_DIRECTION": self.bmap_obj.orbit_direction}
            auto_flip_direction(meta, ax=self.ax_img, print_msg=False)
            self.ax_img.set_xlabel("Range")
            self.ax_img.set_ylabel("Azimuth")
            if self.ax_slide_coh is not None:
                self.sl_last_val = self.sl_coh.val
                self.ax_slide_coh.remove()
                self.ax_slide_coh = None

        par = None
        v_range = None
        cb_ttl = ""
        cmap = None
        if self.rb_par.value_selected == "Velocity":  # show velocity
            v_range = np.max(np.abs(self.vel * self.scale))
            par = self.vel * self.scale
            cb_ttl = f"[{self.vel_scale}/\nyear]"
            cmap = cmc.cm.cmaps["roma"]
        elif self.rb_par.value_selected == "DEM correction":  # show demerr
            v_range = np.max(np.abs(self.demerr))
            par = self.demerr
            cb_ttl = "[m]"
            cmap = cmc.cm.cmaps["vanimo"]

        if self.rb_par.value_selected != "None":
            """idea: show only the points which have any information during this period. Only remove points which have
            no information during the selected time span"""
            # get the indices of the selected time range
            start_idx = self.date_slider.val[0]
            end_idx = self.date_slider.val[1]
            mask_coherent_period = ~np.isnan(self.point_obj.phase)
            mask_coherent_points = mask_coherent_period[:, start_idx:end_idx + 1].astype(np.int8).sum(axis=1) != 0

            self.sc = self.ax_img.scatter(self.point_obj.coord_xy[mask_coherent_points, 1],
                                          self.point_obj.coord_xy[mask_coherent_points, 0],
                                          c=par[mask_coherent_points],
                                          s=5,
                                          cmap=cmap,
                                          vmin=-v_range,
                                          vmax=v_range)
            if self.rb_ccs.value_selected == "CCS":
                # mask and mark all CCS pixels
                self.ax_img.scatter(self.point_obj.coord_xy[self.mask_ccs, 1],
                                    self.point_obj.coord_xy[self.mask_ccs, 0],
                                    c=par[self.mask_ccs],
                                    cmap=cmap,
                                    s=30,
                                    marker='$\u25A1$',
                                    vmin=-v_range,
                                    vmax=v_range)

        self.cb.ax.set_title(cb_ttl, fontsize=self.font_size)
        self.cb = self.fig1.colorbar(self.sc, cax=self.ax_cb, ax=self.ax_img, pad=0.03, shrink=0.8, aspect=10,
                                     orientation='vertical')

        # add back location of selected time series point and current reference
        if self.ts_refpoint_idx is not None:  # initial value is None
            y, x = self.point_obj.coord_xy[self.ts_refpoint_idx, :]
            self.ts_refpoint_marker = self.ax_img.scatter(x, y, marker='^', facecolors='none', edgecolors='k')

        y, x = self.point_obj.coord_xy[self.ts_point_idx, :]
        self.ts_point_marker = self.ax_img.scatter(x, y, facecolors='none', edgecolors='k')

        if not flag_initial_plot:
            self.ax_img.set_xlim(ax_img_xlim)
            self.ax_img.set_ylim(ax_img_ylim)

        plt.draw()

    def updateButtonStatus(self, val: object):  # val seems to be unused, but its necessary for the function to work.
        """Set to true."""
        if self.set_reference_point:
            self.set_reference_point = False
            self.button_mask.color = '1'
        else:
            self.set_reference_point = True
            self.button_mask.color = '0.5'

    def onClick(self, event):
        """Event function to get y/x from button press."""
        if event.inaxes is None:
            return

        if not plt.fignum_exists(self.fig2.number):
            self.initFigureTimeseries()
            plt.show()

        if event.button is MouseButton.RIGHT:
            if event.inaxes == self.ax_img:
                y, x = int(event.ydata + 0.5), int(event.xdata + 0.5)
                idx = self.tree.query([y, x])[-1]
                y, x = self.point_obj.coord_xy[idx, :]

                if self.set_reference_point:  # update reference point
                    if ~self.mask_ccs[idx]:
                        self.logger.warning(msg="The selected point is not used as reference, as it is not continuously"
                                                " coherent. Closest continuously coherent scatterer is used instead.")
                        ccs_idx = self.tree_ccs.query([y, x])[-1]  # index among the CCS points
                        ccs_coord = self.point_obj.coord_xy[self.mask_ccs, :]
                        # query again to get the correct index among all points
                        self.ts_refpoint_idx = self.tree.query(ccs_coord[ccs_idx, :])[-1]
                    else:
                        self.ts_refpoint_idx = idx
                    self.updateReference()
                    self.updateButtonStatus(val=None)
                    # if self.ts_refpoint_marker is not None:  # initial value is None
                    #     self.ts_refpoint_marker.remove()
                    # self.ts_refpoint_marker = self.ax_img.scatter(x, y, marker='^', facecolors='none', edgecolors='k')
                else:
                    self.ts_point_idx = idx

                    self.updateNeighbourhood(val=None)
                    if self.ts_point_marker is not None:  # initial value is None
                        self.ts_point_marker.remove()
                        y, x = self.point_obj.coord_xy[self.ts_point_idx, :]
                    self.ts_point_marker = self.ax_img.scatter(x, y, facecolors='none', edgecolors='k')
                self.plotPointTimeseries(val=None)
        return

    def updateNeighbourhood(self, val):
        """Update the neighbourhood of the selected point."""
        self.neighb_idx, self.neighb_mask = selectNeighbourhood(
            searchtree=self.tree_coord_map,
            coord_map=self.point_obj.coord_map,
            idx=self.ts_point_idx,
            radius=float(self.txt_radius.text)
        )
        self.neighb_markers = plotNeighbourhoodMap(
            ax=self.ax_img,
            neighb_markers=self.neighb_markers,
            neighb_coord_xy=self.point_obj.coord_xy[self.neighb_mask, :]
        )

    def updateReference(self):
        """Change the phase of all points according to the new reference point.

        Update the plot of the velocity and time series.
        """
        self.logger.info(msg="changed reference to ID: {}".format(self.point_obj.point_id[self.ts_refpoint_idx]))
        self.point_obj.phase -= self.point_obj.phase[self.ts_refpoint_idx, :]
        vel, demerr, ref_atmo, coherence, omega, v_hat = ut.estimateParameters(obj=self.point_obj, ifg_space=False)
        self.vel = vel
        self.demerr = demerr
        self.ref_atmo = ref_atmo
        self.plotMap(val=None)

    def updateBackgroundTimeseriesPlot(self, val: object):
        """Update the background of the timeseries plot according to selected time span."""
        self.highlight_time_span.remove()
        self.highlight_time_span = self.ax_ts.axvspan(
            xmin=self.times[int(self.date_slider.val[0])],
            xmax=self.times[int(self.date_slider.val[1])],
            color='grey',
            alpha=0.1
        )
        self.fig2.canvas.draw()

    def prepareTimeseries(self, *, point_idx: int):
        """Prepare phase time series for plotting."""
        # transform phase time series into meters
        resulting_ts = self.point_obj.wavelength / (4 * np.pi) * self.point_obj.phase[point_idx, :]
        cbox_status = self.cbox_par.get_status()
        if not cbox_status[0]:  # Displacement
            resulting_ts = resulting_ts - self.point_obj.ifg_net_obj.tbase * self.vel[point_idx]
        if not cbox_status[1]:  # DEM error
            phase_topo = (self.point_obj.ifg_net_obj.pbase / (self.point_obj.slant_range[point_idx] *
                                                              np.sin(self.point_obj.loc_inc[point_idx])) *
                          self.demerr[point_idx])
            resulting_ts = resulting_ts - phase_topo

        # compute trend
        if self.rb_baselines.value_selected == "Temporal baseline":
            line = self.point_obj.ifg_net_obj.tbase * self.vel[point_idx] + self.ref_atmo[point_idx]
        elif self.rb_baselines.value_selected == "Perpendicular baseline":
            line = (self.point_obj.ifg_net_obj.pbase / (self.point_obj.slant_range[point_idx] *
                                                        np.sin(self.point_obj.loc_inc[point_idx])) *
                    self.demerr[point_idx] + self.ref_atmo[point_idx])
        else:
            line = None
        return resulting_ts * self.scale, line * self.scale

    def plotNeighbourhoodTimeseries(self):
        """Plot timeseries neighbouring points."""
        zorder = 0
        style = {
            "c": "grey",
            "marker": ".",
            "markersize": 0.5,
            "linewidth": 0.5,
            "linestyle": ":",
        }
        for line in self.neighb_ts_lines:
            line[0].remove()

        self.neighb_ts_lines = list()
        for idx in self.neighb_idx:
            resulting_ts = self.prepareTimeseries(point_idx=idx)[0]

            if self.rb_baselines.value_selected == "Temporal baseline":
                self.neighb_ts_lines.append(
                    self.ax_ts.plot(self.times, resulting_ts, zorder=zorder, **style))
                valid_idx = np.where(~np.isnan(resulting_ts))[0]
                self.neighb_ts_lines.append(
                    self.ax_ts.plot([self.times[valid_idx[0]], self.times[valid_idx[-1]]],
                                    [resulting_ts[valid_idx[0]], resulting_ts[valid_idx[-1]]],
                                    '.', zorder=zorder, ** {"c": "grey", "markersize": 3}))

            if self.rb_baselines.value_selected == "Perpendicular baseline":
                # self.neighb_ts_lines.append(
                #     self.ax_ts.plot(self.point_obj.ifg_net_obj.pbase, resulting_ts,
                #                     '.', zorder=zorder, **style))
                pass

    def plotPointTimeseries(self, val: object):  # val seems to be unused, but it is necessary for the function to work.
        """Plot_point_timeseries."""
        ts_query_point, line = self.prepareTimeseries(point_idx=self.ts_point_idx)
        self.ax_ts.set_ylabel(f"Displacement [{self.vel_scale}]")

        if self.rb_fit.value_selected == "Linear fit":
            if self.rb_baselines.value_selected == "Temporal baseline":
                self.lin_fit_line.set_data(self.times, line)
            elif self.rb_baselines.value_selected == "Perpendicular baseline":
                self.lin_fit_line.set_data(self.point_obj.ifg_net_obj.pbase, line)
        else:
            self.lin_fit_line.set_data([], [])

        if self.rb_baselines.value_selected == "Temporal baseline":
            self.ts_line.set_data(self.times, ts_query_point)
            self.ts_line.set_zorder(2)
            self.ax_ts.set_xlabel("Time [years]")
            self.ax_ts.set_xlim(self.times[0], self.times[-1])
        if self.rb_baselines.value_selected == "Perpendicular baseline":
            self.ts_line.set_data(self.point_obj.ifg_net_obj.pbase, ts_query_point)
            self.ts_line.set_zorder(2)
            self.ax_ts.set_xlabel("Perpendicular Baseline [m]")
            self.ax_ts.set_xlim(np.min(self.point_obj.ifg_net_obj.pbase), np.max(self.point_obj.ifg_net_obj.pbase))

        self.text_obj_time.remove()
        point_info = "DEM error: {:.0f} m\nVelocity: {:.0f} {:s}/year".format(
            self.demerr[self.ts_point_idx],
            self.vel[self.ts_point_idx] * self.scale,
            self.vel_scale,
        )
        self.text_obj_time = self.ax_info_box.text(0.5, 0.5, point_info, ha='center', va='center')

        # plot time series of neighbourhood
        self.plotNeighbourhoodTimeseries()

        # set y-lim to [-20, 20] mm except if it exceeds these limits
        y_max = max([0.02 * self.scale, np.nanmax(ts_query_point) + 0.005 * self.scale])
        y_min = min([-0.02 * self.scale, np.nanmin(ts_query_point) - 0.005 * self.scale])
        for line in self.neighb_ts_lines:
            y_min = min([np.nanmin(line[0]._y), y_min])
            y_max = max([np.nanmax(line[0]._y), y_max])
        self.ax_ts.set_ylim(y_min, y_max)

        # update figure
        self.fig1.canvas.draw()
        self.fig2.canvas.draw()


def selectNeighbourhood(searchtree: KDTree, coord_map: np.ndarray, idx: int, radius: float):
    """Select points within a certain radius around a point.

    Parameters
    ----------
    searchtree: KDTree
        searchtree for fast nearest neighbour search
    coord_map: np.ndarray
        map coordinates of the points (dim: no. points x 2)
    idx: int
        index of the point around which the neighbourhood is selected
    radius: float
        radius of the neighbourhood in [m]
    """
    neighb_idx = searchtree.query_ball_point(coord_map[idx, :], r=radius)
    neighb_idx.remove(idx)  # remove the query point
    neighb_mask = np.array([True if i in neighb_idx else False for i in range(len(coord_map))])
    return neighb_idx, neighb_mask


def plotNeighbourhoodMap(ax: plt.Axes, neighb_markers: Optional[PathCollection], neighb_coord_xy: np.ndarray):
    """Plot selected neighbourhood pixels on map."""
    if neighb_markers is not None:
        neighb_markers.remove()
    neighb_markers = ax.scatter(neighb_coord_xy[:, 1], neighb_coord_xy[:, 0], facecolors='none', edgecolors='w')
    return neighb_markers
