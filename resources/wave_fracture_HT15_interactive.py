"""wave_fracture_interactive.py

This Python script provides an interactive plot illustrating wave-ice fracture
as coded into NEMO-SI3, assuming a Bretschneider wave spectrum input. It plots
the wave energy spectrum, the sea surface height (SSH) field computed over the
1D sub-gridscale domain, and the 'fracture histogram' used to update the floe
size distribution. All parameters can be adjusted interactively.

Versions of the software and required packages this has been tested with:

    Python     | 3.12.0
    numpy      | 1.26.0
    matplotlib | 3.8.1


Plots
-----
(a) Wave-energy (Bretschneider) spectrum

(b) Sea surface height corresponding to the wave spectrum in (a). By default
    it is generated with all frequency components having the same phase. The
    points plotted along the SSH indicate where ice breaks assuming it exists
    along the entire 1D subdomain.

(c) The 'fracture histogram': the normalised distribution of all distances
    between the break points indicated in (b), converted to radii of
    resulting floes. Histogram bins are the default floe size distribution
    bins in SI3.


Adjustable parameters
---------------------

Significant wave height, Peak frequency
    These adjust the wave spectrum which is defined here as a Bretschneider
    spectrum.

Mean ice thickness
    Breaking points are calculated using a flexural strain failure criterion,
    where the strain is proportional to the ice thickness. So, adjusting this
    changes the calculation of strain and hence where the ice breaks.

Critical strain
    The threshold strain beyond which ice breaks.

xmax, dx
    The size and resolution of the 1D sub-gridscale domain used to generate
    the realisation of the sea surface height

rmin
    Radius of the smallest floe size assumed to be affected by wave breakup.
    Note this restricts the scale of wavelengths that affects wave breakup,
    but it does not mean that fracture sizes cannot be smaller than 2*rmin.

nf, f0, kf
    Parameters specifying the discretisation of the wave spectrum, as is
    implemented in NEMO: the number of frequences nf, the lowest frequency
    f0, and the constant-ratio kf between spectral classes so that they
    are exponentially spaced.


Button options
--------------
Random phase button
    Re-generate the SSH using a random phase for each frequency component and
    then re-calculate wave fracture. Subsequently adjusting another parameter
    keeps the same set of random phases until this button (or the fixed phase
    button) is clicked again.

Fixed phase button
    Re-generate the SSH using a constant phase (pi) for each frequency
    component and then re-calculate wave fracture.

Monochromatic wave button
    Generate a monochromatic wave instead of using the Bretschneider formula.
    The amplitude and frequency can be adjusted using the significant wave
    height and peak frequency sliders. The random phase button can also be
    used but it makes no difference to the calculation of fracture in this
    case.

    This never enters into NEMO-SI3 but is a useful sanity check/test case
    (e.g., there should be one fracture radius equal to the wavelength,
    the value of which is printed above panel b, divided by four, when
    the strain threshold is exceeded).

Reset button
    Reset to the starting plots.


Other options
-------------
One positional command-line argument may be passed: a color that is passed
to various matplotlib plotting commands, to override the default (tab:red)
used for the break points on panel (b) and fracture histogram in (c).

"""

from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Button, Slider


# Settings for parameter sliders (do not change order):
# |-----------------------------------------------------------------------------------------|
# | param. |    min  | init | max | step  | label                                           |
# |-----------------------------------------------------------------------------------------|
slider_opts = {
    "hs"   : [ [0.00 , 0.30 , 1.00, 0.05 ], r"Significant wave height, $H_\mathrm{s}$ (m)" ],
    "fp"   : [ [0.01 , 0.08 , 0.20, 0.01 ], r"Peak frequency, $f_\mathrm{p}$ (Hz)"         ],
    "ec"   : [ [0.00 , 1.00 , 3.00, 0.25 ], r"Critical strain, $\varepsilon/\varepsilon_0$"],
    "hi"   : [ [0.00 , 0.30 , 1.00, 0.05 ], r"Mean ice thickness, $h_\mathrm{i}$ (m)"      ],
    "xmax" : [ [0.50 , 5.00 , 10.0, 0.50 ], r"$x_\mathrm{max}$ (km)"                       ],
    "dx"   : [ [0.50 , 1.00 , 10.0, 0.50 ], r"$\Delta{}x$ (m)"                             ],
    "rmin" : [ [0.00 , 10.0 , 30.0, 1.00 ], r"$r_\mathrm{min}$ (m)"                        ],
    "nf"   : [ [4    , 25   , 50  , 1    ], r"$n_f$"                                       ],
    "f0"   : [ [0.005, 0.04 , 0.10, 0.005], r"$f_0$ (Hz)"                                  ],
    "kf"   : [ [1.001, 1.10 , 1.20, 0.001], r"$k_f$"                                       ]
}
# |-----------------------------------------------------------------------------------------|

# For constant phase, use this value for all frequencies:
phase_init = np.pi

# Critical strain (value when ec slider set to 1):
ec_0 = 3.e-5

# Bin edges for histogram (hard-coded to match current default FSD bins with nn_nfsd = 12):
hist_bins = np.array([  0.1  ,   2.15987,  10.3142,  28.3481,   59.6495,  107.088, 172.929, 258.78,
                      365.584, 493.639  , 642.659 , 811.847 , 1000.])

# Other constants:
grav = 9.80665

# Plot configuration:
plot_line_color     = "0.1"
plot_marker_size    = 3.
plot_accent_color   = "tab:red" if len(argv) < 2 else argv[1]
button_active_color = "yellow"

# Texts for the random phase (rph), fixed phase (fph), and reset buttons:
fph_button_text = r"$\phi(f)=\pi$"
rph_button_text = r"$\phi(f)\sim\mathrm{U}(0,2\pi)$"
res_button_text = r"$\mathbf{Reset}$"
mon_button_text = u"\u00bd" + r"$H_\mathrm{s}\cos\left(k_\mathrm{p}x+\phi_\mathrm{p}\right)$"

# --------------------------------------------------------------------------- #
# Functions for calculations
# =========================================================================== #

def wfreq(nf_in, f0_in, kf_in):
    """Calculate frequencies and frequency bin widths as calculated in module
    sbcwave.F90: exponentially spaced with a common ratio kf_in, with minimum
    frequency f0_in, and total number of spectral bands nf_in.
    """
    wfreq_out  = np.zeros(nf_in)
    wdfreq_out = np.zeros(nf_in)

    wfreq_out [0] = f0_in * np.sqrt(kf_in)
    wdfreq_out[0] = f0_in * (kf_in - 1.)

    for i in range(1, nf_in):
        wfreq_out[i]  = kf_in * wfreq_out[i-1]
        wdfreq_out[i] = wfreq_out[i] * (np.sqrt(kf_in) - 1./np.sqrt(kf_in))

    return wfreq_out, wdfreq_out


def wlam(wfreq_in):
    """Calculate wavelength(s) from a given input frequency/frequencies. This
    is the dispersion relation for surface deep water gravity waves, as assumed
    in module icewav.F90.
    """
    return grav / (2. * np.pi * wfreq_in**2)


def wspec_bret(hs_in, fp_in, wfreq_in):
    """Calculate the Bretschneider spectrum (as in icewav.F90 function
    wav_spec_bret) from a given significant wave height hs_in, peak
    frequency fp_in, as a function of frequency, wfreq_in.
    """
    return .3125 * hs_in**2 * (fp_in**4 / wfreq_in**5) * np.exp( -1.25 * (fp_in/wfreq_in)**4 )


def x1d(xmax_in, dx_in):
    """Calculate the coordinates of the 1D sub-grid domain for wave fracture,
    from a given domain size (xmax_in) and resolution/spacing (dx_in).
    """
    return np.arange(0., xmax_in + .5*dx_in, dx_in)


def ssh_mono(x1d_in, hs_in, fp_in, fphase_in):
    """Calculate the sea surface height (SSH) along a specified 1D sub-grid domain
    x1d_in, for a monochromatic wave field of frequency fp_in, phase fphase_in,
    and significant wave height hs_in.
    """
    return .5 * hs_in * np.cos(fphase_in + 4. * np.pi**2 * fp_in**2 * x1d_in / grav)


def ssh(x1d_in, wfreq_in, wdfreq_in, wspec_in, fphase_in):
    """Calculate the sea surface height (SSH) along a specified 1D sub-grid domain
    x1d_in, from a wave spectrum wspec_in defined on frequency bands wfreq_in
    with widths wdfreq_in, and phases fphase_in.
    """
    ssh_out = np.zeros(len(x1d_in))
    for i in range(len(wfreq_in)):
        ssh_out[:] += np.sqrt(2. * wspec_in[i] * wdfreq_in[i]) \
                      * np.cos(fphase_in[i] + 4. * np.pi**2 * wfreq_in[i]**2 * x1d_in[:] / grav)
    return ssh_out


def calc_break_locs(x1d_in, ssh_in, dx_in, rmin_in, hi_in, ec_in):
    """Calculate the locations along the 1D sub-grid domain x1d_in, for which
    the sea surface height (SSH) ssh_in is defined, at which sea ice of mean
    thickness hi_in would break according to a flexural strain failure
    criterion. The critical strain threshold is ec_in.

    This function also requires dx_in, the domain spacing, and rmin_in, the
    minimum floe size (as radius) affected by wave fracture. These parameters
    define the size of the 'moving window' used to locate extrema in SSH,
    effectively filtering out waves that have too short a wavelength to have
    any effect on sea ice below the threshold size rmin_in.

    Calculations follow the same method implemented in SI3, subroutine
    wav_frac_dist of module icewav.F90, where possible, except some of the
    indexing differs due to language differences and the window size is
    approximated from rmin_in and dx_in (in SI3 rmin is specified as an
    integer number of dx), so there may be some minor round-off error.
    """

    nx = len(x1d_in)

    # Locate the extrema in SSH; create boolean arrays to track whether each point
    # along x1d_in corresponds to an minima or to a maxima:
    llmax = np.zeros(nx, dtype=bool)
    llmin = np.zeros(nx, dtype=bool)

    # Extrema are defined as maxima or minima over a 'window' of size corresponding
    # to floes of diameter = 2*rmin.
    j_win_half = int(rmin_in // dx_in)  # half of the window size, corresponds to
                                        # nn_ice_wav_rmin in icewav.F90 code

    for jx in range(j_win_half, nx - j_win_half):

        ixlo = jx - j_win_half
        ixhi = jx + j_win_half + 1

        llmin[jx] = np.argmin(ssh_in[ixlo:ixhi]) == j_win_half
        llmax[jx] = np.argmax(ssh_in[ixlo:ixhi]) == j_win_half

    llext = np.logical_or(llmax, llmin)

    # Loop over extrema points and determine whether ice breaks there or not.
    # Create a boolean array, the output of this function, to track which
    # points break:
    breaks_out = np.zeros(nx, dtype=bool)

    for jx in range(j_win_half + 1, nx - j_win_half - 1):
        if llext[jx]:

            # Current point is an extremum. Locate the nearest extrema on
            # each side, indices ixlo and ixhi. Then, if this 'triplet' of
            # extrema alternates (max., min., max.) or (min., max., min.),
            # then calculate the strain at the centre point (jx).

            ixlo = -1  # 'missing' value if no triplet can be found
            ixhi = -1

            # Find nearest extremum on the left:
            for jy in range(jx-1, -1, -1):
                if llext[jy]:
                    ixlo = jy
                    break

            # Find nearest extremum on the right:
            for jy in range(jx+1, nx, 1):
                if llext[jy]:
                    ixhi = jy
                    break

            if ixlo >= 0 and ixhi >= 0:
                if (   (llmin[ixlo] and llmax[jx] and llmin[ixhi])
                    or (llmax[ixlo] and llmin[jx] and llmax[ixhi])):

                    # Calculate strain across this triplet of points
                    #
                    # There is a finite difference approximation for the second
                    # second derivative of ssh_in, with respect to x1d_in, that
                    # simplifies to the expression calculated below.
                    #
                    dxlo = x1d_in[jx]   - x1d_in[ixlo]
                    dx   = x1d_in[ixhi] - x1d_in[ixlo]
                    dxhi = x1d_in[ixhi] - x1d_in[jx]

                    e = abs(.5 * hi_in * (ssh_in[ixhi]*dxlo - ssh_in[jx]*dx + ssh_in[ixlo]*dxhi)
                                       / (dxlo*dx*dxhi) )

                    # Save whether this point breaks or not:
                    breaks_out[jx] = e >= ec_in

    return breaks_out


def calc_frac_sizes(x1d_in, breaks_in):
    """Calculate an array of fracture lengths from a given 1D sub-grid domain,
    x1d_in, and the corresponding boolean array for which locations break
    (breaks_in) as calculated in function calc_break_locs().
    """
    return x1d_in[breaks_in][1:] - x1d_in[breaks_in][:-1]


def calc_frac_histogram(frac_sizes_in):
    """Calculate the normalised 'fracture histogram' from a given array of
    fracture sizes, frac_sizes_in, calculated in function calc_frac_sizes().
    The histogram bins are currently hardcoded (and set as a global array
    hist_bins).
    """
    if len(frac_sizes_in) > 0:
        return np.histogram(frac_sizes_in, density=True, bins=hist_bins)[0]
    else:
        return np.zeros(len(hist_bins)-1)

# --------------------------------------------------------------------------- #

if __name__ == "__main__":

    print(__doc__)

    # ----------------------------------------------------------------------- #
    # Calculate initial data
    # ======================================================================= #

    # Calculate all required fields from initial values of sliders; these will
    # be used as global variables in the functions to update plots via widgets:

    data_wfreq, data_wdfreq = wfreq(slider_opts["nf"][0][1], slider_opts["f0"][0][1],
                                    slider_opts["kf"][0][1])

    data_is_mono = False
    data_wspec   = wspec_bret(slider_opts["hs"][0][1] , slider_opts["fp"][0][1], data_wfreq)
    data_x1d     = x1d(slider_opts["xmax"][0][1]*1000., slider_opts["dx"][0][1])
    data_x1d_km  = data_x1d * 1.e-3

    # Create phases for maximum number of frequencies selectable, so even if parameter nf
    # (number of frequencies) changes we still do not need to calculate random phases again
    # (initially, phases are fixed):
    nf_max      = slider_opts["nf"][0][2]
    data_phase  = phase_init * np.ones(nf_max)
    data_ssh    = ssh(data_x1d, data_wfreq, data_wdfreq, data_wspec, data_phase)

    data_breaks = calc_break_locs(data_x1d, data_ssh,
                                  rmin_in=slider_opts["rmin"][0][1],
                                  dx_in  =slider_opts["dx"  ][0][1],
                                  hi_in  =slider_opts["hi"  ][0][1],
                                  ec_in  =slider_opts["ec"  ][0][1] * ec_0)

    data_frac_sizes = calc_frac_sizes(data_x1d, data_breaks)
    data_frac_hist  = calc_frac_histogram(data_frac_sizes)


    # ----------------------------------------------------------------------- #
    # Set up figure with initial plots
    # ======================================================================= #

    plt.rcdefaults()

    fig = plt.figure(figsize=(9.6, 4.8))
    fig.canvas.manager.set_window_title("Wave Fracture")
    axd = fig.subplot_mosaic("AB;CB", width_ratios=[1.,2.5])

    line_wspec, = axd["A"].plot(data_wfreq, data_wspec, color=plot_line_color,
                                marker="o", ms=plot_marker_size, label=r"$S$")

    line_fp     = axd["A"].axvline([slider_opts["fp"][0][1]], linestyle=":",
                                   color=plot_line_color, label=r"$f_\mathrm{p}$")

    line_ssh,   = axd["B"].plot(data_x1d_km, data_ssh, color=plot_line_color,
                                label="SSH", zorder=1)

    sctr_break  = axd["B"].scatter(data_x1d_km[data_breaks], data_ssh[data_breaks],
                                   color=plot_accent_color, s=plot_marker_size**2,
                                   label="Break points", zorder=2)

    hist_bin_widths = hist_bins[1:] - hist_bins[:-1]

    bar_frac    = axd["C"].bar(hist_bins[:-1], data_frac_hist, width=hist_bin_widths,
                               align="edge", facecolor=plot_accent_color, edgecolor="w")

    axd["A"].set_xlim(xmin=0.)
    axd["A"].set_ylim(0., .15)
    axd["A"].set_xlabel(r"$f$ (Hz)")
    axd["A"].set_ylabel(r"$S(H_\mathrm{s},f_\mathrm{p})$ (m$^2$ Hz$^{-1}$)")
    axd["A"].set_title("(a) Wave spectrum", fontweight="bold")
    axd["A"].legend(loc="upper right")

    axd["B"].set_xlim(0., slider_opts["xmax"][0][2])
    axd["B"].set_ylim(-.6, .6)
    axd["B"].set_xlabel(r"$x$ (km)")
    axd["B"].set_ylabel(r"SSH (m)")
    axd["B"].set_title("(b) Sea surface height", fontweight="bold")
    axd["B"].legend(loc="upper right", ncols=2, columnspacing=.4, handletextpad=.25, scatterpoints=3)

    ssh_text = axd["B"].annotate(r"$\lambda_\mathrm{p}=%.1f$ m" % wlam(slider_opts["fp"][0][1]),
                                 (.01, 1.02), ha="left", va="bottom", xycoords="axes fraction")

    axd["C"].set_xlim(0., 500.)
    axd["C"].set_ylim(0., .03)
    axd["C"].set_xlabel(r"Fracture size (m)")
    axd["C"].set_ylabel(r"Frequency density (m$^{-1}$)")
    axd["C"].set_title("(c) Fracture histogram", fontweight="bold")

    for k in axd.keys():
        axd[k].set_label(axd[k].get_title())

    # Create a function to return a string for the histogram statistics, so that
    # it can be used to update the text on the plot in the widget functions:
    def stat_str(x):
        """Returns a string with various statistical values of an input array x
        to be used on plots (intended for x = fracture lengths).
        """
        if len(x) > 0:
            xstr = (  f"Range: [{np.min(x):.1f}, {np.max(x):.1f}] m" + "\n"
                    + f"Median = {np.median(x):.1f} m" + "\n"
                    + f"IQR = {np.subtract(*np.percentile(x, [75, 25])):.1f} m")
            if len(x) > 1:
                xstr += "\n" + f"Skewness = {np.mean(((x - np.mean(x)) / np.std(x))**3):.1f}"
        else:
            xstr = ""

        return xstr

    hist_stat_text = axd["C"].annotate(stat_str(data_frac_sizes), (.96, .95),
                                       ha="right", va="top", xycoords="axes fraction",
                                       fontsize="smaller")

    # Set overall layout before adjusting below to make space for and create widgets:
    fig.tight_layout()

    # ----------------------------------------------------------------------- #
    # Create slider widgets
    # ======================================================================= #

    # Make space for widgets by adjusting the SSH axis, then set positions of
    # the slider widgets. This part is somewhat trial-and-error to get the
    # right placement, but it works well enough.
    #
    # Set the vertical displacement of the bottom edge of the SSH axes (dy_axB).
    #
    # Then draw a background box above which all the sliders and buttons sit;
    # this will be positioned to line up with the right hand edge of the SSH
    # axes, a specified y-coordinate of the bottom edge (box_y0), a height set
    # to match the lower edge of the SSH axes but with a manually-set margin
    # (dy_box), and similarly for the width (adjusted with dx_box).
    #
    # All of these measures are in figure coordinates.
    #
    box_y0 = .03    # y-coordinate of box lower edge
    dy_axB = .35    # vertical displacement of SSH plot axes
    dx_box = -.05   # horizontal adjustment of box left edge
    dy_box = -.125  # vertical adjustment of box top edge

    axd["B"].set_position([axd["B"].get_position().x0, axd["B"].get_position().y0 + dy_axB,
                           axd["B"].get_position().width,
                           axd["B"].get_position().height - dy_axB])

    # Fix positions of each plot (so they don't automatically move if adjusting
    # limits etc. interactively, although it doesn't seem to stop this problem
    # if using the interactive 'configure subplots', alas):
    for k in axd.keys():
        axd[k].set_position(axd[k].get_position(), which="both")

    # Add the background box first. The widgets are positioned relative to this box,
    # so calculate and save the coordinates, width, and height:
    box_xy0 = (axd["B"].get_position().x0 + dx_box, box_y0)
    box_wid = axd["B"].get_position().x1 - box_xy0[0]
    box_hei = axd["B"].get_position().y0 + dy_box - box_xy0[1]

    fig.patches.extend([mpatches.Rectangle(box_xy0, width=box_wid, height=box_hei,
                                           facecolor=".975", edgecolor="k", zorder=-10,
                                           linewidth=.8, transform=fig.transFigure)])

    # Coordinates of sliders: they are arranged as a row of four, then a leftward-
    # displaced array of 3 rows by 2 columns (to make space for buttons).
    #
    # Use coordinates as fraction of box width of height (i.e., similar to figure
    # coordinates but scaled to the background box rather than whole figure canvas).
    #
    # Specify the vertical margin and height of each slider, from which their
    # vertical coordinates can be calculated:
    #
    smg = .06  # y-margin between slider/background box
    sh0 = .09  # height of each slider

    sy0 = (1. - smg - sh0)  # y-coordinate, first/top row (lower edge of slider)
    sdy = (sy0 - smg) / 6.  # y-spacing

    # Specify the x-coordinates (again, as fraction of background box width)
    #
    # Slider widths: one value sw[0] for the top row of four, then sw[1] for
    # the first column of the bottom three rows, and then sw[2] for the
    # second column:
    #
    sw = [.52, .14, .14]

    # Slider positions (edges of sliders themselves, not their labels):
    sx = [.4, .5, .4 + sw[0] - sw[2]]

    # Keyword arguments passed to all Slider instances:
    slider_kw = {"initcolor": "1", "track_color": ".8", "facecolor": ".2"}

    def _create_slider(j, k):
        """Helper function to create sliders (don't want to save these in a dictionary like
        above for the axes as it is slower to repeatedly access a dictionary interactively).
        This takes an index j (numbering slider axes from 0-3 for top four rows, 4-6 for
        first column of smaller sliders, 7-9 for second column of smaller slides) and a
        parameter key (k, as used for keys in dictionary in slider_opts) as inputs and
        returns a Slider instance positioned appropriately on the figure.
        """

        jj = 0 if j < 4 else (1 if j < 7 else 2)  # index of sw and sx

        ax_j = fig.add_axes([box_xy0[0] + sx[jj] * box_wid,
                             box_xy0[1] + (sy0 - j*sdy + 3.*sdy*(j>6)) * box_hei,
                             sw[jj]*box_wid, sh0*box_hei])

        ax_j.set_label(f"<Slider widget: {k}>")

        return Slider(ax=ax_j, label=slider_opts[k][1], valmin=slider_opts[k][0][0],
                      valmax=slider_opts[k][0][2], valinit=slider_opts[k][0][1],
                      valstep=slider_opts[k][0][3], **slider_kw)

    slider_hs   = _create_slider(0, "hs")
    slider_fp   = _create_slider(1, "fp")
    slider_ec   = _create_slider(2, "ec")
    slider_hi   = _create_slider(3, "hi")
    slider_xmax = _create_slider(4, "xmax")
    slider_dx   = _create_slider(5, "dx")
    slider_rmin = _create_slider(6, "rmin")
    slider_nf   = _create_slider(7, "nf")
    slider_f0   = _create_slider(8, "f0")
    slider_kf   = _create_slider(9, "kf")

    # Set axes for the buttons (in the 'recessed' area before the bottom 3 rows of sliders).
    # Buttons are arranged in a row of two by two. Use a similar system for sliders: x- and
    # y-coordinates as fraction of background box width/height:
    #
    bx0 = .014  # 1st column (left edge of button)
    bx1 = .140  # 2nd column (left edge of button)
    by0 = .26   # 1st row (bottom edge of button)
    by1 = .06   # 2nd row (bottom edge of button)
    bw0 = .11   # 1st column, width
    bw1 = .22   # 2nd column, width
    bh0 = .14   # button heights (same for both rows)

    ax_fph = fig.add_axes([box_xy0[0] + bx0*box_wid, box_xy0[1] + by0*box_hei, bw0*box_wid, bh0*box_hei])
    ax_rph = fig.add_axes([box_xy0[0] + bx1*box_wid, box_xy0[1] + by0*box_hei, bw1*box_wid, bh0*box_hei])
    ax_res = fig.add_axes([box_xy0[0] + bx0*box_wid, box_xy0[1] + by1*box_hei, bw0*box_wid, bh0*box_hei])
    ax_mon = fig.add_axes([box_xy0[0] + bx1*box_wid, box_xy0[1] + by1*box_hei, bw1*box_wid, bh0*box_hei])

    ax_fph.set_label("<Button widget: fixed phase>")
    ax_rph.set_label("<Button widget: random phase>")
    ax_res.set_label("<Button widget: reset>")
    ax_mon.set_label("<Button widget: mono wave field>")

    # Keyword arguments passed to all buttons and colors for showing whether the random
    # or fixed phase button is active or inactive:
    button_kw             = {"hovercolor": ".975"}
    button_inactive_color = "white"

    # Create button instances:
    button_rph = Button(ax_rph, rph_button_text, color=button_inactive_color, **button_kw)
    button_fph = Button(ax_fph, fph_button_text, color=button_active_color  , **button_kw)
    button_res = Button(ax_res, res_button_text, color=button_inactive_color, **button_kw)
    button_mon = Button(ax_mon, mon_button_text, color=button_inactive_color, **button_kw)


    # ----------------------------------------------------------------------- #
    # Create 'action' functions to be called when sliders are adjusted or
    # when buttons are pressed:
    # ======================================================================= #

    # To minimise re-calculations, there are three functions for the sliders:
    #
    #     (1) update_spectrum : hs  , fp, nf, f0, kf
    #     (2) update_ssheight : xmax, dx
    #     (3) update_fracture : ec  , hi, rmin
    #
    # The parameters listed above for each function are those for which their
    # sliders directly call when adjusted. Function (1) calls (2), and function
    # (2) calls (3). This way only the necessary re-calculations are done for
    # each parameter (e.g., ec does not change the wave spectrum of SSH, so
    # only the data relating to fracture calculations need to be re-done,
    # whereas changing hs changes the wave spectrum, which also changes the SSH
    # and hence also the fracture.

    def update_spectrum(val, phase_random=False, phase_reset=False):
        """Actions for updating the wave spectrum plot. This function is used directly
        by the hs, fp, nf, f0, and kf sliders, as changing those parameters directly
        affects the wave spectrum. It is not used at all by the remaining parameter sliders.
        """

        global data_phase, nf_max, data_wfreq, data_wdfreq, data_wspec

        # Set the phase to either random or same initial fixed value, or leave unchanged:
        if phase_random:
            data_phase = 2. * np.pi * np.random.random(nf_max)
        elif phase_reset:
            data_phase = phase_init * np.ones(nf_max)

        # If both phase_random and phase_reset are false, leave phases unchanged
        # so that if a random phase was previously selected, the same set of
        # random phases is kept when adjusting another parameter)

        data_wfreq, data_wdfreq = wfreq(slider_nf.val, slider_f0.val, slider_kf.val)

        if data_is_mono:
            data_wspec = np.zeros(slider_nf.val)
            line_fp.set_linestyle("-")
        else:
            data_wspec = wspec_bret(slider_hs.val, slider_fp.val, data_wfreq)
            line_fp.set_linestyle("--")

        line_wspec.set_data(data_wfreq, data_wspec)
        line_fp.set_xdata([slider_fp.val])

        update_ssheight(val)


    def update_ssheight(val):
        """Actions for updating the sea surface height plot. This function is used directly
        by the xmax and dx sliders, as changing those parameters does not affect the wave
        spectrum. It is not used at all by the sliders for ec, hi, and rmin, and for all
        remaining parameters this function is used indirectly (via the update_spectrum
        function).
        """

        global data_x1d   , data_x1d_km, data_is_mono, data_wfreq, \
               data_dwfreq, data_wspec , data_phase  , data_ssh

        data_x1d    = x1d(slider_xmax.val*1000., slider_dx.val)
        data_x1d_km = data_x1d * 1.e-3

        if data_is_mono:
            data_ssh = ssh_mono(data_x1d, slider_hs.val, slider_fp.val, data_phase[0])
        else:
            data_ssh = ssh(data_x1d, data_wfreq, data_wdfreq, data_wspec, data_phase)

        ssh_text.set_text(r"$\lambda_\mathrm{p}=%.1f$ m" % wlam(slider_fp.val))
        line_ssh.set_data(data_x1d_km, data_ssh)

        update_fracture(val)


    def update_fracture(val):
        """Actions for updating the fracture points and histogram. This function is used
        directly by the ec, hi, and rmin sliders, as this is the only part affected by
        those parameters. It is used indirectly (via other functions) for the sliders for
        all remaining parameters.
        """

        global data_x1d, data_x1d_km, data_ssh, data_breaks, data_frac_sizes, data_frac_hist

        data_breaks = calc_break_locs(data_x1d, data_ssh, rmin_in=slider_rmin.val,
                                      dx_in=slider_dx.val, hi_in=slider_hi.val,
                                      ec_in=slider_ec.val * ec_0)

        data_frac_sizes = calc_frac_sizes    (data_x1d, data_breaks)
        data_frac_hist  = calc_frac_histogram(data_frac_sizes)

        sctr_break.set_offsets(np.array([data_x1d_km[data_breaks], data_ssh[data_breaks]]).T)

        [bar.set_height(data_frac_hist[i]) for i, bar in enumerate(bar_frac)]
        hist_stat_text.set_text(stat_str(data_frac_sizes))

        # This function is always the last step (all parameter changes result in
        # a change to the fracture distribution), so refresh the figure here:
        fig.canvas.draw_idle()


    def random_phase(event):
        """Actions for when 'random phase' button is clicked."""
        button_rph.color = button_active_color
        button_fph.color = button_inactive_color
        update_spectrum(event, phase_random=True)


    def fixed_phase(event):
        """Actions for when 'fixed phase' button is clicked."""
        button_fph.color = button_active_color
        button_rph.color = button_inactive_color
        update_spectrum(event, phase_reset=True)


    def toggle_mono_wave_field(event):
        """Actions for when the 'mono wave field' button is clicked."""
        global data_is_mono
        data_is_mono = not data_is_mono
        button_mon.color = button_active_color if data_is_mono else button_inactive_color
        update_spectrum(event)


    def reset(val):
        """Actions for when 'reset' button is clicked (reset everything to initial
        values, including resetting the phase and deactivating monochromatic wave).
        """

        button_fph.color = button_active_color
        button_rph.color = button_inactive_color
        button_mon.color = button_inactive_color

        global data_is_mono
        data_is_mono = False

        slider_hs.reset()
        slider_fp.reset()
        slider_ec.reset()
        slider_hi.reset()
        slider_xmax.reset()
        slider_dx.reset()
        slider_rmin.reset()
        slider_nf.reset()
        slider_f0.reset()
        slider_kf.reset()

        update_spectrum(val, phase_reset=True)


    # Associate sliders and buttons with the appropriate action functions:
    slider_hs.on_changed(update_spectrum)
    slider_fp.on_changed(update_spectrum)
    slider_ec.on_changed(update_fracture)
    slider_hi.on_changed(update_fracture)
    slider_xmax.on_changed(update_ssheight)
    slider_dx.on_changed(update_ssheight)
    slider_rmin.on_changed(update_fracture)
    slider_nf.on_changed(update_spectrum)
    slider_f0.on_changed(update_spectrum)
    slider_kf.on_changed(update_spectrum)

    button_fph.on_clicked(fixed_phase)
    button_rph.on_clicked(random_phase)
    button_res.on_clicked(reset)
    button_mon.on_clicked(toggle_mono_wave_field)

    # Show the figure interactively:
    plt.show()

