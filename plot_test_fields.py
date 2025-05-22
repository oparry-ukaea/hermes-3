from boututils.datafile import DataFile
from matplotlib import pyplot as plt


def plot_field(var, data, ax, x_index=0):
    ax.imshow(data[x_index, :, :])
    ax.set_title(f"2D slice of {var}")
    ax.set_xlabel("y")
    ax.set_ylabel("z")


def plot_test_fields(f):
    with DataFile(f) as d:
        nk = len(d.keys())
        if nk > 0:
            print(f"Plotting {nk} fields from {f}")
            fig, ax = plt.subplots(nk)
            for ii, var in enumerate(d.keys()):
                plot_field(var, d[var], ax[ii])
            fig.savefig("/home/oparry/code/hermes-3/tmp.png")
        else:
            print(f"No fields in {f}")


plot_test_fields("/home/oparry/code/hermes-3/tests/unit/reactions/Drec.nc")
