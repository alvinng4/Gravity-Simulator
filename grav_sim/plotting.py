"""
Plotting module for the grav_sim package.
"""

from pathlib import Path
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from . import utils

SOLAR_SYSTEM_COLORS = {
    "Sun": "orange",
    "Mercury": "slategrey",
    "Venus": "wheat",
    "Earth": "skyblue",
    "Mars": "red",
    "Jupiter": "darkgoldenrod",
    "Saturn": "gold",
    "Uranus": "paleturquoise",
    "Neptune": "blue",
    "Moon": "grey",
    "Pluto": None,
    "Ceres": None,
    "Vesta": None,
}


def set_3d_axes_equal(ax: plt.Axes) -> None:
    """
    Make axes of 3D plot have equal scale

    Parameters
    ----------
    ax : matplotlib axis
        The axis to set equal scale

    Reference
    ---------
    karlo, https://stackoverflow.com/questions/13685386/how-to-set-the-equal-aspect-ratio-for-all-axes-x-y-z
    """

    x_limits = ax.get_xlim3d()  # type: ignore
    y_limits = ax.get_ylim3d()  # type: ignore
    z_limits = ax.get_zlim3d()  # type: ignore

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5 * max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])  # type: ignore
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])  # type: ignore
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])  # type: ignore


def plot_2d_trajectory(
    sol_state: np.ndarray,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$",
    ylabel: str = "$y$",
    title: Optional[str] = None,
    marker: str = "o",
    markersize: int = 6,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot 2D trajectory of objects

    Parameters
    ----------
    sol_state : np.ndarray
        Solution state of the system
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect="equal")

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    for i in range(sol_state.shape[1]):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None

        traj = ax.plot(
            sol_state[:, i, 1],
            sol_state[:, i, 2],
            color=colors_i,
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i, 1],
            sol_state[-1, i, 2],
            color=traj[0].get_color(),
            label=labels_i,
            marker=marker,
            markersize=markersize,
        )

    if title is not None:
        ax.set_title(title)

    if legend:
        fig.legend(loc="center right", borderaxespad=0.2)
        fig.tight_layout()

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")


def plot_3d_trajectory(
    sol_state: np.ndarray,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    xlabel: str = "$x$",
    ylabel: str = "$y$",
    zlabel: str = "$z$",
    title: Optional[str] = None,
    marker: str = "o",
    markersize: int = 6,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot 3D trajectory of objects

    Parameters
    ----------
    sol_state : np.ndarray
        Solution state of the system
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    xlabel : str, optional
        Label of x-axis, by default "$x$ (AU)"
    ylabel : str, optional
        Label of y-axis, by default "$y$ (AU)"
    zlabel : str, optional
        Label of z-axis, by default "$z$ (AU)"
    marker : str, optional
        Marker for the last position, by default "o"
    markersize : int, optional
        Marker size for the last position, by default 6
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([1.0, 1.0, 1.0])  # type: ignore
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)  # type: ignore

    for i in range(sol_state.shape[1] // 6):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None

        traj = ax.plot(
            sol_state[:, i * 3],
            sol_state[:, i * 3 + 1],
            sol_state[:, i * 3 + 2],
            color=colors_i,
        )
        # Plot the last position with marker
        ax.plot(
            sol_state[-1, i * 3],
            sol_state[-1, i * 3 + 1],
            sol_state[-1, i * 3 + 2],
            color=traj[0].get_color(),
            label=labels_i,
            marker=marker,
            markersize=markersize,
        )

    set_3d_axes_equal(ax)

    if title is not None:
        ax.set_title(title)

    if legend:
        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
        fig.subplots_adjust(
            right=0.7
        )  # Adjust the right boundary of the plot to make room for the legend

    fig.tight_layout()

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")

def plot_quantity_against_time(
    quantity: np.ndarray,
    sol_time: np.ndarray,
    is_log_y: bool = False,
    title: Optional[str] = None,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    colors: Optional[list[str]] = None,
    labels: Optional[list[str]] = None,
    legend: bool = False,
    save_fig: bool = False,
    save_fig_path: Optional[str | Path] = None,
) -> None:
    """Plot a quantity against time

    Parameters
    ----------
    quantity : np.ndarray
        Quantity to be plotted
    sol_time : np.ndarray
        Solution time
    is_log_y : bool, optional
        Flag to check whether y-axis is in log scale, by default False
    title : Optional[str], optional
        Title of the plot, by default None
    xlabel : Optional[str], optional
        Label of x-axis, by default None
    ylabel : Optional[str], optional
        Label of y-axis, by default None
    colors : Optional[list[str]], optional
        Colors of the trajectories, by default None
    labels : Optional[list[str]], optional
        Labels of the objects, used for legend, by default None
    legend : bool, optional
        Flag to check whether to show legend, by default False
    save_fig : bool, optional
        Flag to check whether to save the figure, by default False
    save_fig_path : Optional[str], optional
        Path to save the figure, by default None
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)

    if len(quantity.shape) == 1:
        quantity = quantity.reshape(-1, 1)

    for i in range(quantity.shape[1]):
        if colors is not None:
            colors_i = colors[i]
        else:
            colors_i = None

        if labels is not None:
            labels_i = labels[i]
        else:
            labels_i = None
    
        ax.plot(sol_time, quantity[:, i], color=colors_i, label=labels_i)

    if is_log_y:
        ax.set_yscale("log")
    if title is not None:
        ax.set_title(title)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    if legend:
        ax.legend(loc="center right", bbox_to_anchor=(1.325, 0.5))
        fig.subplots_adjust(
            right=0.7
        )  # Adjust the right boundary of the plot to make room for the legend

    if save_fig:
        plt.savefig(save_fig_path, dpi=300)
    else:
        plt.show()
    plt.close("all")
