#!/usr/bin/env python3

import xhermes

# Following two lines needed so all variables are shown when printing the Dataset
import xarray as xr

xr.set_options(display_max_rows=1000)
# Set better figure size
from matplotlib import pyplot as plt

plt.rcParams["figure.figsize"] = (24, 8)

import os.path


def get_run_dir(eg_name):
    repo_root = os.path.dirname(os.path.dirname(__file__))
    return os.path.join(repo_root, "example-runs", eg_name)


def blob2dteti_anim(rd, **kwargs):
    ds = xhermes.open(rd, **kwargs)
    print(ds)

    # Get rid of size-1 parallel direction
    ds = ds.squeeze()

    # Make an animation
    # Note: saving a gif can be slow. Comment out `save_as` argument to disable.
    ds.bout.animate_list(
        ["Ne", "Vort", "phi"],
        ncols=3,
        show=True,
        save_as="blob2d-teti-anim",
    )


def blob2dteti_plot(d, **kwargs):
    from boutdata import collect
    import matplotlib.pyplot as plt

    n = collect("Ne", path=rd)

    fig, axs = plt.subplots(1, 3, figsize=(9, 3))

    for ax, time in zip(axs, [15, 30, 45]):
        ax.contourf(n[time, :, 0, :].T, 50)

    plt.tight_layout()

    plt.savefig("blob2d-te-ti.png")

    plt.show()


# blob2dteti_anim(rd, unnormalise=False)
blob2dteti_plot(get_run_dir("blob2d-te-ti"), unnormalise=False)
