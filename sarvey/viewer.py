#!/usr/bin/env python

# SARvey - A multitemporal InSAR time series tool for the derivation of displacements.
#
# Copyright (C) 2021-2024 Andreas Piter (IPI Hannover, piter@ipi.uni-hannover.de)
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
from typing import Any
from logging import Logger
import matplotlib.cm as cm
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from matplotlib import colormaps, widgets
from matplotlib.backend_bases import MouseButton
from matplotlib.colors import Normalize
import numpy as np
from scipy.spatial import KDTree
import datetime

from mintpy.objects.colors import ColormapExt
from mintpy.utils import readfile
from mintpy.utils.plot import auto_flip_direction

from sarvey.objects import AmplitudeImage, Points, BaseStack
import sarvey.utils as ut


def plotIfgs(*, phase: np.ndarray, coord: np.ndarray, spatial_ref_idx: int = None, ttl: str = None, cmap: str = "cmy"):
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
        colormap, use "cmy" for wrapped phase data (default) or "?" for unwrapped or residual phase
    """
    if cmap == "cmy":
        cmap = ColormapExt('cmy').colormap
    else:
        cmap = plt.get_cmap(cmap)

    num_ifgs = phase.shape[1]
    min_val = np.min(phase)
    max_val = np.max(phase)
    fig, axs = plt.subplots(np.ceil(np.sqrt(num_ifgs + 1)).astype(np.int32),
                            np.ceil(np.sqrt(num_ifgs + 1)).astype(np.int32))
    sc = None
    for i, ax in enumerate(axs.flat):
        if i < num_ifgs:
            sc = ax.scatter(coord[:, 1], coord[:, 0], c=phase[:, i],
                            vmin=min_val, vmax=max_val, s=1, cmap=cmap)
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
                unit: str = None, s: float = 5.0, cmap: colormaps = colormaps["jet_r"], symmetric: bool = False,
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
        colormap (default: "jet_r")
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
        sc = ax.scatter(coord[:, 1], coord[:, 0], c=value, s=s, cmap=plt.get_cmap(cmap),
                        vmin=-v_range, vmax=v_range)
    else:
        sc = ax.scatter(coord[:, 1], coord[:, 0], c=value, s=s, cmap=plt.get_cmap(cmap), **kwargs)
    cb = plt.colorbar(sc, ax=ax, pad=0.03, shrink=0.5)
    cb.ax.set_title(unit)
    ax.set_title(ttl)
    plt.tight_layout()
    return fig, ax, cb


def plotColoredPointNetwork(*, x: np.ndarray, y: np.ndarray, arcs: np.ndarray, val: np.ndarray, ax: plt.Axes = None,
                            linewidth: float = 2, cmap_name: str = "seismic", clim: tuple = None):
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
    cmap_name: str
        name of the colormap (default: "seismic")
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
        fig = plt.figure(figsize=[15, 5])
        ax = fig.add_subplot()
    else:
        fig = ax.get_figure()
    ax.scatter(x, y, s=3.5, c=np.ones_like(x))

    if clim is None:
        norm = Normalize(vmin=min(val), vmax=max(val))
    else:
        norm = Normalize(vmin=clim[0], vmax=clim[1])

    mapper = cm.ScalarMappable(norm=norm, cmap=cm.get_cmap(cmap_name))
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
        self.plotPointTimeseries(val=None)  # just any point
        self.fig1.canvas.mpl_connect('button_press_event', self.onClick)
        plt.show()

    def initFigureMap(self):
        """InitFigureMap."""
        self.fig1 = plt.figure()
        self.ax_img = self.fig1.subplots(1, 1)

        self.ax_cb = self.fig1.add_axes([0.93, 0.6, 0.015, 0.15])  # (left, bottom, width, height)
        self.cb = self.fig1.colorbar(self.sc,
                                     cax=self.ax_cb,
                                     ax=self.ax_img,
                                     pad=0.03,
                                     shrink=0.8,
                                     aspect=10,
                                     orientation='vertical')

        # add button to select reference point
        self.set_reference_point = False
        self.ax_button = self.fig1.add_axes([0.125, 0.9, 0.1, 0.08])  # (left, bottom, width, height)
        self.button_mask = widgets.Button(ax=self.ax_button, label='Select\nReference', image=None, color='1')
        self.button_mask.on_clicked(self.updateButtonStatus)

        # add radiobutton to select parameter
        self.ax_radio_par = self.fig1.add_axes([0.225, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.rb_par = widgets.RadioButtons(self.ax_radio_par, labels=['Velocity', 'DEM error', 'None'], active=0)
        self.rb_par.on_clicked(self.plotMap)

        # add radiobutton to select background image
        self.ax_radio_backgr = self.fig1.add_axes([0.425, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.rb_backgr = widgets.RadioButtons(self.ax_radio_backgr, labels=['Amplitude', 'DEM', 'Coherence', 'None'],
                                              active=0)
        self.rb_backgr.on_clicked(self.plotMap)

        # add info box with info about velocity and DEM error of selected pixel
        self.ax_info_box = self.fig1.add_axes([0.625, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.text_obj_time = self.ax_info_box.text(0.1, 0.1, "")
        self.ax_info_box.set_xticks([], [])
        self.ax_info_box.set_yticks([], [])

        # add variable for axis of slider controlling the visualized coherence background image
        self.ax_slide_coh = None
        self.sl_last_val = 0.0
        self.sl_coh = None

    def initFigureTimeseries(self):
        """InitFigureTimeseries."""
        self.fig2 = plt.figure(figsize=(15, 5))
        self.ax_ts = self.fig2.subplots(1, 1)

        # add radiobutton for fitting linear model
        self.ax_radio_fit = self.fig2.add_axes([0.125, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.rb_fit = widgets.RadioButtons(self.ax_radio_fit, labels=['None', 'Linear fit'], active=0)

        # add radiobutton for selecting baseline type
        self.ax_radio_baselines = self.fig2.add_axes([0.325, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.rb_baselines = widgets.RadioButtons(
            self.ax_radio_baselines,
            labels=['Temporal baseline', 'Perpendicular baseline'],
            active=0
        )

        # add check box for removing phase due to parameters
        self.ax_cbox_par = self.fig2.add_axes([0.525, 0.9, 0.2, 0.08])  # (left, bottom, width, height)
        self.cbox_par = widgets.CheckButtons(
            self.ax_cbox_par,
            ["Velocity", "DEM error"],
            actives=[True, False]
        )
        self.rb_fit.on_clicked(self.plotPointTimeseries)
        self.rb_baselines.on_clicked(self.plotPointTimeseries)
        self.cbox_par.on_clicked(self.plotPointTimeseries)

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
                self.ax_slide_coh = self.fig1.add_axes([0.425, 0.85, 0.2, 0.03])  # (left, bottom, width, height)
                self.sl_coh = widgets.Slider(self.ax_slide_coh,
                                             label='Coherence',
                                             valmin=0.0,
                                             valmax=1.0,
                                             valinit=self.sl_last_val,
                                             valfmt="%.1f")

            self.ax_img.imshow(self.temp_coh_img,
                               cmap=cm.grayC,
                               vmin=np.round(self.sl_coh.val, decimals=1),
                               vmax=1)
            meta = {"ORBIT_DIRECTION": self.bmap_obj.orbit_direction}
            auto_flip_direction(meta, ax=self.ax_img, print_msg=False)
            self.ax_img.set_xlabel("Range")
            self.ax_img.set_ylabel("Azimuth")
        if self.rb_backgr.value_selected == "None":
            self.ax_img.imshow(np.ones_like(self.height, dtype=np.int8), cmap=plt.cm.get_cmap("gray"), vmin=0, vmax=1)
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
        if self.rb_par.value_selected == "Velocity":  # show velocity
            v_range = np.max(np.abs(self.vel * self.scale))
            par = self.vel * self.scale
            cb_ttl = f"[{self.vel_scale}/\nyear]"
        elif self.rb_par.value_selected == "DEM error":  # show demerr
            v_range = np.max(np.abs(self.demerr))
            par = self.demerr
            cb_ttl = "[m]"

        if self.rb_par.value_selected != "None":
            self.sc = self.ax_img.scatter(self.point_obj.coord_xy[:, 1],
                                          self.point_obj.coord_xy[:, 0],
                                          c=par,
                                          s=5,
                                          cmap=cm.roma,
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
                    self.ts_refpoint_idx = idx
                    self.updateReference()
                    self.updateButtonStatus(val=None)
                    # if self.ts_refpoint_marker is not None:  # initial value is None
                    #     self.ts_refpoint_marker.remove()
                    # self.ts_refpoint_marker = self.ax_img.scatter(x, y, marker='^', facecolors='none', edgecolors='k')
                else:
                    self.ts_point_idx = idx

                    if self.ts_point_marker is not None:  # initial value is None
                        self.ts_point_marker.remove()
                        y, x = self.point_obj.coord_xy[self.ts_point_idx, :]
                    self.ts_point_marker = self.ax_img.scatter(x, y, facecolors='none', edgecolors='k')
                self.plotPointTimeseries(val=None)
        return

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

    def plotPointTimeseries(self, val: object):  # val seems to be unused, but its necessary for the function to work.
        """Plot_point_timeseries."""
        self.ax_ts.cla()

        # transform phase time series into meters
        resulting_ts = self.point_obj.wavelength / (4 * np.pi) * self.point_obj.phase[self.ts_point_idx, :]
        cbox_status = self.cbox_par.get_status()
        if not cbox_status[0]:  # Displacement
            resulting_ts = resulting_ts - self.point_obj.ifg_net_obj.tbase * self.vel[self.ts_point_idx]
        if not cbox_status[1]:  # DEM error
            phase_topo = (self.point_obj.ifg_net_obj.pbase / (self.point_obj.slant_range[self.ts_point_idx] *
                                                              np.sin(self.point_obj.loc_inc[self.ts_point_idx])) *
                          self.demerr[self.ts_point_idx])
            resulting_ts = resulting_ts - phase_topo

        self.ax_ts.set_ylabel(f"Displacement [{self.vel_scale}]")

        # add trend
        if self.rb_fit.value_selected == "Linear fit":
            if self.rb_baselines.value_selected == "Temporal baseline":
                line = self.point_obj.ifg_net_obj.tbase * self.vel[self.ts_point_idx] + self.ref_atmo[self.ts_point_idx]
                self.ax_ts.plot(self.times, line * self.scale, '-k')
            elif self.rb_baselines.value_selected == "Perpendicular baseline":
                line = (self.point_obj.ifg_net_obj.pbase / (self.point_obj.slant_range[self.ts_point_idx] *
                                                            np.sin(self.point_obj.loc_inc[self.ts_point_idx])) *
                        self.demerr[self.ts_point_idx] + self.ref_atmo[self.ts_point_idx])
                self.ax_ts.plot(self.point_obj.ifg_net_obj.pbase, line * self.scale, '-k')

        # set y-lim to [-20, 20] mm except if it exceeds this scale
        y_max = max([0.02, resulting_ts.max() + 0.005])
        y_min = min([-0.02, resulting_ts.min() - 0.005])

        self.ax_ts.set_ylim(y_min * self.scale, y_max * self.scale)
        if self.rb_baselines.value_selected == "Temporal baseline":
            self.ax_ts.plot(self.times, resulting_ts * self.scale, '.')
            self.ax_ts.set_xlabel("Time [years]")
        if self.rb_baselines.value_selected == "Perpendicular baseline":
            self.ax_ts.plot(self.point_obj.ifg_net_obj.pbase, resulting_ts * self.scale, '.')
            self.ax_ts.set_xlabel("Perpendicular Baseline [m]")

        self.text_obj_time.remove()
        point_info = "DEM error: {:.0f} m\nVelocity: {:.0f} {:s}/year".format(
            self.demerr[self.ts_point_idx],
            self.vel[self.ts_point_idx] * self.scale,
            self.vel_scale,
        )
        self.text_obj_time = self.ax_info_box.text(0.5, 0.5, point_info, ha='center', va='center')

        # update figure
        self.fig1.canvas.draw()
        self.fig2.canvas.draw()
