"""fsd_cat_bounds.py

Originally used to prototype and test the category definition user interface in
the NEMO-SI3 implementation of the prognostic FSD, this script plots the FSD
categories for the various options that can calculate them internally in SI3.

Four options are available: user-defined in the namelist (not shown here),
uniform spacing, Gaussian spacing, and exponential spacing. The latter three
are all determined by specifying a minimum and maximum floe size to be resolved,
and a number of floe size categories. For the Gaussian and exponential spacing,
an additional parameter, called 's' here, is available to adjust the shape to make
it more linear (i.e., less refined spacing at smaller categories).

Note that increasing s beyond 1 for the exponential case quickly leads to extremely
small categories that might not even be able to be resolved. In SI3 (icefsd.F90,
subroutine fsd_initbounds), a warning is raised whenever the spacing becomes smaller
than 1cm (in any initialisation case) with suggestions to (and how to) adjust.

Increasing s > 1 in the Gaussian case does not make much difference to the shape
[there is a limit -- see docstring of function bounds_gaussian()].

These methods allow online production of a variety of floe size discretisations,
and can approximate well enough the hard-coded limits used in Icepack (optionally
plotted here) and both 'partitions' used in the model of Zhang et al. (2015, JGR:O).
Close approximations to these are given by:

    Icepack bounds, n = 12                  --> Gaussian    with s ~ 10
    Icepack bounds, n = 24                  --> Exponential with s ~ 0.75
    Zhang et al. (2015) partition 1, n = 12 --> Gaussian    with s ~ 10

(the second partition of Zhang et al. 2015 has uniform spacing, which is trivially
reproduced and not included in this script).

This script also creates a figure showing the floe welding array (called 'floe_iweld'
in SI3 and just 'iweld' in Icepack) corresponding to the category that each pair of
floe-floe welding interactions results in.

Use --help for usage/options.

"""

from argparse import ArgumentParser

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate


# Plotting parameters (using keys 'hc' for hard-coded limits, 'un' for uniform,
#                      'gs' for Gaussian, and 'ex' for exponential)
colors = {"hc": "tab:grey", "un": "tab:blue", "gs": "tab:green", "ex": "tab:orange"}


def bounds_uniform(n, rmin, rmax):
    """Returns n+1 limits for n floe size categories that are uniformly spaced
    covering the range rmin to rmax inclusive.
    """
    return np.linspace(rmin, rmax, n+1)


def bounds_gaussian(n, rmin, rmax, s=1.):
    """Returns n+1 limits for n floe size categories with a Gaussian spacing and
    spanning the range rmin to rmax inclusive. Specifically, the limits L(i) are
    formulated similarly to how the ice thickness category limits are defined by
    Hibler (1980; Mon. Weather Rev.; Appendix C), as cited by Zhang et al. (2015):

                ( rmin                                        i = 0
        L(i) = <
                ( L(i-1) + k * [ 1 - exp( -(i/(s*n))^2 ) ]    i = 1..n

    where k is a constant that is calculated to ensure L(n) = rmax and s is a shape
    parameter that can be thought of as the standard deviation and controls the
    degree of non-linearity in the limits. The exponent includes a factor of n
    multiplying s, so that changing the number of categories without changing the
    overall range (rmin and/or rmax) keeps the same overall shape (note that the
    expression does not depend on n for i = 0 and i = n).

    Here, s can be varied down to a limit of 0 (but not exactly), which makes the
    spacing uniform. The default of 1 makes it close to as 'fully curved' as it
    can be; increasing it further does not make much difference as there is a
    limit to how 'curved' the Gaussian can be with a constraint on L(n) = rmax
    (mathematically, increasing s decreases k and the two changes compensate in
    the limit of large s).
    """

    s = max(1.e-10, s)    # avoid division by zero

    # Determine k by using L(n) = rmax = rmin + k * sum(1 - exp(...))
    # from the recursive formula for L(i) and then rearranging for k:
    k = 0.
    for i in range(0, n):
        k += 1. - np.exp(-((n-i)/(s*n))**2)
    k = (rmax - rmin) / k

    lims = np.zeros(n+1)
    lims[0] = rmin
    for i in range(1, n+1):
        lims[i] = lims[i-1] + k * (1. - np.exp(-(i/(s*n))**2))

    return lims


def bounds_exponential(n, rmin, rmax, s):
    """Returns n+1 limits for n floe size categories with exponentially-increasing
    spacing and spanning the range rmin to rmax inclusive. The limits L(i) are
    defined to have a similar form to the Gaussian profile:

                ( rmin                            i = 0
        L(i) = <
                ( L(i-1) + k * exp( 10*s*i/n )    i = 1..n

    where k is a constant that is calculated to ensure L(n) = rmax and s is a shape
    parameter that controls the degree of non-linearity in the limits. The exponent
    includes a factor of n dividing s, so that changing the number of categories without
    changing the overall range (rmin and/or rmax) keeps the same overall shape (noting
    that the expression does not depend on n for i = 0 and i = n). There is also an
    ad-hoc factor of 10 to ensure that with default parameters, the floe spacing at small
    categories is not too small, i.e., setting the default scale to not be 'too exponential'.

    Here, s can be varied down to a limit of 0, which makes the spacing uniform. It can
    be increased to anything, but this can easily result in spacings that are too small
    to be resolved. There is not really anything that can be done to avoid it so a warning
    is added in SI3 if the spacing ever becomes smaller than 1 cm.
    """

    # Determine k by using L(n) = rmax = rmin + k * exp(...)
    # from the recursive formula for L(i) and then rearranging for k:
    k = 0.
    for i in range(1, n+1):
        k += np.exp(10.*s*i/n)
    k = (rmax - rmin) / k

    lims = np.zeros(n+1)
    lims[0] = rmin
    for i in range(1, n+1):
        lims[i] = lims[i-1] + k * np.exp( 10.*s*i/n )

    return lims


def calc_iweld(rl, rc, ru):
    """Calculate floe welding array corresponding to 'floe_iweld' in SI3 and representing
    the category that each pair of category-category interactions results in due to
    welding. Here, floe area is taken as the floe radius squared and the 'shape'
    parameter is ignored (or implicitly set to 1, but it doesn't make a difference anyway).

    The welding array is determined as it is currently done in SI3/Icepack by considering
    the sum of areas of floes at the centres of categories j1 and j2, and finding the
    category j3 = iweld(j1,j2) whose limits contain that net floe size.
    """

    nfsd = len(rl)
    iweld = np.zeros((nfsd,nfsd), dtype=int)

    for j1 in range(nfsd):
        for j2 in range(nfsd):
            aweld = rc[j1]**2 + rc[j2]**2  # ~ area of welded floe
            for j3 in range(nfsd-1):  # find welded category
                if aweld >= rl[j3]**2 and aweld < ru[j3]**2:
                    iweld[j1,j2] = j3
            if aweld >= rl[nfsd-1]**2:  # check for largest category (no upper limit)
                iweld[j1,j2] = nfsd-1

    iweld += 1  # so indices are Fortran style / correspond to category number

    return iweld


def main():

    prsr = ArgumentParser(usage="Plot FSD category limits")
    prsr.add_argument("-n"          , type=int  , default=12  , help="Number of categories")
    prsr.add_argument("-s"          , type=float, default=1.  , help="Non-linearity parameter")
    prsr.add_argument("-m", "--rmin", type=float, default=None, help="Minimum floe size")
    prsr.add_argument("-M", "--rmax", type=float, default=None, help="Maximum floe size")
    prsr.add_argument("--hard-coded", type=str  , default="ip", help="Hard-coded limits to plot",
                                      choices=["ip", "z15", "none"])
    prsr.add_argument("-w", "--weld", type=str  , default="hc", help="Which limits for welding array",
                                      choices=["hc", "un", "gs", "ex", "none"])

    prsr.add_argument("-p", "--print-lims", action="store_true",
                      help="Print limit values in a table")
    prsr.add_argument("-z", "--z15-lims", action="store_true",
                      help="Compare with Zhang et al. (2015) limits, not Icepack")

    cmd = prsr.parse_args()

    lims = {}  # store limits in dict with keys 'hc', 'un', 'gs', and 'ex'

    # Set hard-coded limits for comparison
    if cmd.hard_coded == "none":
        n = cmd.n
        rmin = 0.1   if cmd.rmin is None else cmd.rmin
        rmax = 1000. if cmd.rmax is None else cmd.rmax
        lims["hc"] = np.nan * np.zeros(n+1)
    else:
        if cmd.hard_coded == "z15":
            # Limits for 'partition 1' (Gaussian spaced) in Zhang et al. (2015):
            lims["hc"] = np.array([-0.1, 0.1, 10.2, 40.2, 99.8, 199.1, 347.9, 556.,
                                   833.1, 1189.2, 1633.9, 2176.8, 2827.7])
        else:
            # Hard-coded limits in Icepack:
            lims["hc"] = np.array([.0665000, 5.3103085, 14.2865861, 29.0576686, 52.4122136,
                                   87.8691405, 139.5184700, 211.6357520, 308.0372740,
                                   431.2030590, 581.2772250, 755.1410470, 945.8128340,
                                   1343.5444600, 1822.6536400, 2472.6136100, 3354.3498800,
                                   4550.5141300, 6173.2316400, 8374.6117000, 11361.0059000,
                                   15412.3510000, 20908.4095000, 28364.3675000, 38479.1270000])

        n    = min(cmd.n, len(lims["hc"])-1)
        lims["hc"] = lims["hc"][:n+1]
        rmin = lims["hc"][0]  if cmd.rmin is None else cmd.rmin
        rmax = lims["hc"][-1] if cmd.rmax is None else cmd.rmax

    # X-axes for plotting:
    cat    = np.arange(1., n+1, 1.)
    catlim = np.arange(.5, n+1, 1.)

    # Save lower limit (rl), centre of category (rc), and upper limit (ru):
    rl = {}  ;  rc = {}  ;  ru = {}

    # ---- Hard-coded bounds and radii arrays:
    rl["hc"] = lims["hc"][:n]
    ru["hc"] = lims["hc"][1:n+1]
    rc["hc"] = .5 * (rl["hc"] + ru["hc"])

    # ---- Bounds and radii arrays calculated from formula:
    lims["un"] = bounds_uniform(n, rmin, rmax)
    rl["un"]   = lims["un"][:n]
    ru["un"]   = lims["un"][1:n+1]
    rc["un"]   = .5 * (rl["un"] + ru["un"])

    lims["gs"] = bounds_gaussian(n, rmin, rmax, cmd.s)
    rl["gs"]   = lims["gs"][:n]
    ru["gs"]   = lims["gs"][1:n+1]
    rc["gs"]   = .5 * (rl["gs"] + ru["gs"])

    lims["ex"] = bounds_exponential(n, rmin, rmax, cmd.s)
    rl["ex"]   = lims["ex"][:n]
    ru["ex"]   = lims["ex"][1:n+1]
    rc["ex"]   = .5 * (rl["ex"] + ru["ex"])

    if cmd.print_lims:
        if cmd.hard_coded == "none":
            headers = ["i", "Gaussian", "Exponential", "Uniform"]
            rows = []
            for j in range(n+1):
                rows.append([j] + [lims[x][j] for x in ["gs", "ex", "un"]])
        else:
            headers = ["i", "Icepack" if cmd.hard_coded == "ip" else "Zhang et al. (2015)",
                       "Gaussian", "Exponential", "Uniform"]
            rows = []
            for j in range(n+1):
                rows.append([j] + [lims[x][j] for x in ["hc", "gs", "ex", "un"]])

        print("\n" + ("-"*44))
        print(  f"Category limits computed for n = {n}\nr_min = {100.*rmin:.1f} cm, "
              + f"r_max = {rmax/1000.:.3f} km, s = {cmd.s:.2f}")
        print(("="*44) + "\n")
        print(tabulate(rows, headers=headers) + "\n" + ("-")*44)

    # ----------------------------------------------------------------------- #
    # Plot the limits computed from the various functions (and hard-coded)
    # ======================================================================= #
    mpl.rcParams["font.size"] = 14.

    # ax.scatter() keyword arguments for plotting centres rc and limits rl/ru:
    centre_kw = {"marker": "o", "s": 10. , "clip_on": False, "edgecolor": "none"}
    limits_kw = {"marker": "+", "s": 100., "clip_on": False}

    fig1, ax1 = plt.subplots(figsize=(6.4,5.25))
    fig1.canvas.manager.set_window_title("Floe size category definitions")

    for k, label in zip(["un"      , "gs"       , "ex"],
                        ["Uniform" , "Gaussian" , "Exponential"]):

        ax1.plot(   catlim, lims[k]/1000., color=colors[k], label=label)
        ax1.scatter(cat-.5, rl[k]  /1000., color=colors[k], **limits_kw)
        ax1.scatter(cat   , rc[k]  /1000., color=colors[k], **centre_kw)
        ax1.scatter(cat+.5, ru[k]  /1000., color=colors[k], **limits_kw)

    if cmd.hard_coded != "none":

        ax1.plot(catlim, lims["hc"]/1000., color=colors["hc"], linestyle="--",
                 label="Zhang et al. (2015)" if cmd.hard_coded == "z15"
                       else "Icepack (hardcoded)")

        ax1.scatter(cat-.5, rl["hc"]/1000., color=colors["hc"], **limits_kw)
        ax1.scatter(cat   , rc["hc"]/1000., color=colors["hc"], **centre_kw)
        ax1.scatter(cat+.5, ru["hc"]/1000., color=colors["hc"], **limits_kw)

    # For legend:
    ax1.scatter(np.nan, np.nan, color="k", label="Category limits" , **limits_kw)
    ax1.scatter(np.nan, np.nan, color="k", label="Category centres", **centre_kw)

    # Format plot:
    ax1.set_title("Floe size category definitions")
    ax1.set_xlabel("Category")
    ax1.set_ylabel("Radius (km)")
    ax1.set_xlim(.5, n+.5)
    ax1.set_ylim(ymin=0)
    ax1.legend(fontsize="x-small")

    ax1.xaxis.set_major_locator(mpl.ticker.FixedLocator(cat))
    ax1.xaxis.set_minor_locator(mpl.ticker.FixedLocator(catlim))
    ax1.tick_params(axis="x", which="major", size=0)
    ax1.tick_params(axis="x", which="minor", size=5)

    if n > 18:
        ax1.tick_params(axis="x", which="major", labelsize="x-small")

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_position(("outward", 10))
    ax1.spines["bottom"].set_position(("outward", 10))

    # Set figure layout and add a text label with the parameters at the bottom:
    fig1.tight_layout()
    fig1.subplots_adjust(bottom=.2, top=.925)

    fig1.text(.025, .025,  (r"$n=%i$  |  $r_\mathrm{min}=%.1f$ cm  |  " % (n, 100.*rmin))
                         + (r"$r_\mathrm{max}=%.3f$ km  |  $s=%.2f$"    % (rmax/1000., cmd.s)),
              ha="left", va="bottom", fontsize="x-small")


    # ----------------------------------------------------------------------- #
    # Plot the floe welding array for the specified method
    # ======================================================================= #
    if cmd.weld != "none" and not (cmd.weld == "hc" and cmd.hard_coded == "none"):
        fig2, ax2 = plt.subplots(figsize=(4.8, 5.25))
        fig2.canvas.manager.set_window_title("Floe welding array")

        # Calculate floe welding array for specified set of limits. Create two
        # copies, one for the pcolormesh where the invalid values are masked out
        # and one retaining the actual array values for the text labels:
        iweld_pcm = calc_iweld(rl[cmd.weld], rc[cmd.weld], ru[cmd.weld])
        iweld_txt = iweld_pcm.copy()

        # Mask out invalid values in iweld_pcm:
        for i in range(n):
            for j in range(n):
                if i > j:
                    iweld_pcm[i,j] = -1  # one side of diagonal
                elif i == j == iweld_pcm[i,j] - 1:
                    iweld_pcm[i,j] = -1  # cats. welding with themselves into same cat.

        # Create a discrete colormap with one color per category index:
        cmap = mpl.colormaps["Spectral"]
        cmap.set_under("lightgrey")  # for masked out values (-1)
        norm = mpl.colors.BoundaryNorm(catlim, 256)

        pcm = ax2.pcolormesh(catlim, catlim, iweld_pcm, edgecolor="w",
                             cmap=cmap, norm=norm)

        # Line plot to mark the 'discrete diagonal':
        ax2.step(catlim, catlim, where="pre", color="k", clip_on=False)

        # Label each cell with the actual array value (even if masked)
        for i in range(n):
            for j in range(n):
                ax2.annotate(f"{iweld_txt[i,j]}", (j+1, i+1),
                             color="k" if iweld_pcm[i,j] > 0 else "tab:grey",
                             ha="center", va="center", fontweight="normal",
                             fontsize="x-small" if n<=16 else 6)

        # Format axes:
        kw = {"aspect": "equal"}
        for x, j in zip("xy", "12"):
            kw[f"{x}lim"]   = (.5, n+.5)
            kw[f"{x}label"] = f"Category {j}"
            getattr(ax2, f"{x}axis").set_major_locator(mpl.ticker.FixedLocator(cat))

        ax2.set(**kw)
        ax2.set_title("Floe welding array")
        ax2.tick_params(which="major", axis="both", length=0,
                        labelsize="small" if n<=16 else 8)

        for spine in ax2.spines:
            ax2.spines[spine].set_visible(False)

        # Set figure layout and add text description at the bottom:
        fig2.tight_layout()
        fig2.subplots_adjust(top=.925, bottom=.25)

        if cmd.weld == "hc":
            y = "Icepack hard-coded" if cmd.hard_coded == "ip" else "Zhang et al. (2015)"
        else:
            y  = "uniformly" if cmd.weld == "un" else (
                 "Gaussian"  if cmd.weld == "gs" else "exponentially")
            y += "-spaced"

        x  =  "Plot shows the category that results from the welding of floes in category 1 "
        x += f"and\ncategory 2 using the {y} limits. Grey cells represent interactions\n"
        x +=  "not computed to avoid double counting (above the 'diagonal'/black line) or if\n"
        x +=  "categories 1 and 2 are the same and weld to within the same category (along\n"
        x += u"'diagonal'). This array corresponds to the variable 'floe_iweld' in SI\u00b3 "
        x +=  "('iweld' in\nIcepack) and is calculated in the same way here."

        fig2.text(0.025, 0.015, x, ha="left", va="bottom", fontsize="xx-small")

    plt.show()


if __name__ == "__main__":
    main()
