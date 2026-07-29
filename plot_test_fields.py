from matplotlib import pyplot as plt
from glob import glob
import numpy as np
import os.path
import netCDF4 as nc

src_fld_names=['density_source','momentum_source','energy_source']
slices=["xy","xz","yz"]
midpoint_indices=dict(xy=3,xz=2,yz=1)

def plot_field(var, data, ax, slice="xy"):

    index = midpoint_indices.get(slice,0)
    dmin = np.min(data)
    eps=1e-30
    offset = (1+1e-6)*(eps-dmin)
    if slice=="xy":
        logdata =np.log10(1+data[:, :, index]) if dmin>0 else np.log10(data[:, :, index]+offset)
    elif slice=="yz":
        logdata =np.log10(1+data[index, :, :]) if dmin>0 else np.log10(1+data[index, :, :]+offset)
    elif slice=="xz":
        logdata =np.log10(1+data[:, index, :]) if dmin>0 else np.log10(1+data[:, index, :]+offset)

    im = ax.imshow(logdata)
    #im = ax.contourf(logdata)
    cbar = plt.colorbar(im)
    if dmin > 0:
        cbar.set_label(f"log10 {var}")
    else:
        cbar.set_label(f"log10 {var}+{offset:.2E}")
        
    # N.B. data shape results in imshow/contourf plotting the axes backwards, so reverse the labels
    ax.set_xlabel(slice[1])
    ax.set_ylabel(slice[0])

def plot_sp_fields(sp_name, sp_data, lbl):
    fld_names = [v for v in sp_data.variables if sp_data[v].size>1]

    print(f"Plotting {sp_name}")
    # Report scalar field vals
    for fn in [ varname for varname in sp_data.variables if sp_data[varname].size==1]:
        print(f"  Ignoring scalar field {fn} = {sp_data[fn].getValue().max()}")

    nf = len(fld_names)

    if nf == 0:
        return

    fig, ax = plt.subplots(nrows=nf,ncols=3,figsize=(12,3*nf+2))
    fig.suptitle(f"Species [{sp_name}]", fontsize=16)
    for ifld,fld_name in enumerate(fld_names):
        for idim in range(3):
            slice=slices[idim]
            fld_data = sp_data[fld_name]
            try:
                plot_field(fld_name, fld_data, ax[ifld][idim],slice=slice)
            except TypeError:
                plot_field(fld_name, fld_data, ax[ifld][idim],slice=slice)
    fig.tight_layout(pad=1.0)
    fig.savefig(f"/home/oparry/code/hermes-3/{lbl}_{sp_name}.png")

def plot_test_fields(f):
    lbl = os.path.splitext(os.path.basename(f))[0]
    with nc.Dataset(f) as d:
        all_sp = d['species']
        print(f"Plotting {len(all_sp.groups)} species from {f}")
        for sp_name,sp_data in all_sp.groups.items():
            plot_sp_fields(sp_name,sp_data,lbl)

data_files = glob("./tests/unit/reactions/*.nc")
for data_file in data_files:
    plot_test_fields(data_file)
