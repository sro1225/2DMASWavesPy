def plot_cmp_section(config):

    import glob
    import pickle
    import os
    import re
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.ndimage import gaussian_filter
    
    plotcfg = config["plotting"]
    MODEL_DIR = (config["output"]["model_dir"])
    site = config["project"]["site"]
    profile = config["project"]["profile"]
    # Input parameter
    DZ = plotcfg["dz"]
    MAX_DEPTH = (plotcfg["max_depth"])
    SIGMA_VERTICAL = (plotcfg["sigma_vertical"])
    SIGMA_HORIZONTAL = (plotcfg["sigma_horizontal"])
    CMAP_VS = plotcfg["cmap_vs"]
    CMAP_STD = plotcfg["cmap_std"]
    ###########################################################################
    # FIND AND LOAD MODELS
    ##########################################################################
    files = sorted(glob.glob(os.path.join(MODEL_DIR,"__CMP_*_mean.pkl")))

    if len(files) == 0:
        raise RuntimeError("No pickle files found.")

    cmp_locations = []
    models = []

    print(f"Loading {len(files)} models...")

    for fname in files:
        match = re.search(r"CMP_([0-9.]+)",os.path.basename(fname))

        if match is None:
            print(f"Skipping {fname}")
            continue

        cmp_loc = float(match.group(1))

        with open(fname, "rb") as f:
            model = pickle.load(f)

            cmp_locations.append(cmp_loc)
            models.append(model)

    cmp_locations = np.asarray(cmp_locations)

    # sort by CMP location
    sort_idx = np.argsort(cmp_locations)

    cmp_locations = cmp_locations[sort_idx]
    models = [models[i] for i in sort_idx]

    ##########################################################################
    # BUILD CELL EDGES
    ##########################################################################

    x_edges = np.zeros(len(cmp_locations) + 1)

    if len(cmp_locations) > 1:

        x_edges[1:-1] = (
            cmp_locations[:-1] + cmp_locations[1:]
        ) / 2.0

        dx_first = cmp_locations[1] - cmp_locations[0]
        dx_last = cmp_locations[-1] - cmp_locations[-2]

        x_edges[0] = cmp_locations[0] - dx_first / 2.0
        x_edges[-1] = cmp_locations[-1] + dx_last / 2.0
    else:
        x_edges[0] = cmp_locations[0] - 0.5
        x_edges[1] = cmp_locations[0] + 0.5

    #####################################################################
    # COMMON DEPTH GRID
    #####################################################################

    #zmax = max(np.max(model["z"]) for model in models)
    #  MAX_DEPTH = 33.0  # m
    zmax = MAX_DEPTH
    z_grid = np.arange(0, zmax + DZ, DZ)
    z_edges = np.empty(len(z_grid) + 1)
    z_edges[1:-1] = (z_grid[:-1] + z_grid[1:]) / 2
    z_edges[0] = z_grid[0] - DZ/2
    z_edges[-1] = z_grid[-1] + DZ/2
    
    nz = len(z_grid)
    nx = len(cmp_locations)
    
    #####################################################################
    # BUILD VELOCITY AND UNCERTAINTY GRIDS
    #####################################################################

    vs_grid = np.full((nz, nx), np.nan)
    std_grid = np.full((nz, nx), np.nan)
    
    for ix, model in enumerate(models):

        beta = np.asarray(model["beta"])
        
        z_interfaces = np.asarray(model["z"])
        
        beta_low = np.asarray(model["beta_low"])
        beta_up  = np.asarray(model["beta_up"])
        
        beta_std = (beta_up - beta_low) / 2.0
        
        # Layer tops and bottoms
        layer_tops = np.r_[0.0, z_interfaces]
        layer_bottoms = np.r_[z_interfaces, zmax]
        
        for vs, vs_std, top, bottom in zip(
                beta,
                beta_std,
                layer_tops,
                layer_bottoms):
            
            mask = (
                (z_grid >= top) &
                (z_grid < bottom)
            )
            
            vs_grid[mask, ix] = vs
            std_grid[mask, ix] = vs_std
        
    ###########################################################################
    # SMOOTHED VERSIONS
    ###########################################################################

    vs_smooth = gaussian_filter(
        vs_grid,
        sigma=(SIGMA_VERTICAL, SIGMA_HORIZONTAL)
    )

    std_smooth = gaussian_filter(
        std_grid,
        sigma=(SIGMA_VERTICAL, SIGMA_HORIZONTAL)
    )

    ##########################################################################
    # COLOUR LIMITS
    ##########################################################################

    vmin_vs = np.nanmin(vs_grid)
    vmax_vs = np.nanmax(vs_grid)

    vmin_std = np.nanmin(std_grid)
    vmax_std = np.nanmax(std_grid)

    ##########################################################################
    # FIGURE
    ##########################################################################

    fig, axes = plt.subplots(
        2, 2,
        figsize=(16, 10),
        sharex=True,
        sharey=True,
        constrained_layout=True
    )

    fig.suptitle(
        f"{site} - {profile}\nPseudo-2D Vs Model and Uncertainty",
        fontsize=18,
        fontweight="bold"
    )


    ###########################################################################
    # BLOCKY VELOCITY
    ###########################################################################
    print("vs_grid:", vs_grid.shape)
    print("vs_smooth:", vs_smooth.shape)

    print("std_grid:", std_grid.shape)
    print("std_smooth:", std_smooth.shape)

    # axis definitions
    for ax in axes.flat:
        ax.invert_yaxis()

    for ax in axes[0,:]:
        ax.tick_params(
            top=True,
            labeltop=True,
            bottom=False,
            labelbottom=False
        )

    pcm1 = axes[0, 0].pcolormesh(
        x_edges,
        z_edges,
        vs_grid,
        shading="auto",
        cmap=CMAP_VS,
        vmin=vmin_vs,
        vmax=vmax_vs,
    )

    axes[0, 0].invert_yaxis()
    axes[0, 0].set_title("Blocky Vs Model")
    axes[0, 0].set_ylabel("Depth (m)")

    # Overlay interface uncertainties on smoothed Vs panel
    ax = axes[0,0]
    for x, model in zip(cmp_locations, models):
        z = np.asarray(model["z"])
        z_low = np.asarray(model["z_low"])
        z_up = np.asarray(model["z_up"])
        ax.errorbar(
            np.full(len(z), x),
            z,
            yerr=[
                z - z_low,
                z_up - z
            ],
            fmt='k.',
            markersize=4,
            alpha=0.8,
            lw=1.5,
            capsize=3,
            capthick=1.5
        )

        ####################################################################
        # SMOOTHED VELOCITY
        ####################################################################

    pcm2 = axes[0, 1].pcolormesh(
        x_edges,
        z_edges,
        vs_smooth,
        shading="auto",
        cmap=CMAP_VS,
        vmin=vmin_vs,
        vmax=vmax_vs,
    )

    axes[0, 1].set_title("Smoothed Vs Model")

    # Overlay interface uncertainties on smoothed Vs panel
    ax = axes[0,1]

    for x, model in zip(cmp_locations, models):
        z = np.asarray(model["z"])
        z_low = np.asarray(model["z_low"])
        z_up = np.asarray(model["z_up"])

        ax.errorbar(
            np.full(len(z), x),
            z,
            yerr=[
                z - z_low,
                z_up - z
            ],
            fmt='k.',
            markersize=4,
            alpha=0.8,
            lw=1.5,
            capsize=3,
            capthick=1.5
        )

    ######
    # Add some contour lines
    # Automatic contour lines
    #####

    #levels = np.arange(100, 1200, 100)
    # exclude half space
    mask_depth = z_grid <= 20

    # thicker black contours underneath
    axes[0,1].contour(
        cmp_locations,
        z_grid[mask_depth],
        vs_smooth[mask_depth, :],
        levels=10,
        colors='black',
        linewidths=1.5,
        alpha=0.7
    )

    # thin white contours on top
    cs = axes[0,1].contour(
        cmp_locations,
        z_grid[mask_depth],
        vs_smooth[mask_depth, :],
        levels=10,
        colors='white',
        linewidths=0.75
    )

    axes[0,1].clabel(
        cs,
        fmt='%d',
        fontsize=8,
        colors='white'
    )
    ###########################################################################
    # BLOCKY UNCERTAINTY
    ###########################################################################
    
    pcm3 = axes[1, 0].pcolormesh(
        x_edges,
        z_edges,
        std_grid,
        shading="auto",
        cmap=CMAP_STD,
        vmin=vmin_std,
        vmax=vmax_std,
    )

    axes[1, 0].set_title("Blocky Vs Uncertainty")
    #axes[1, 0].set_xlabel("Distance along profile (m)")
    axes[1, 0].set_ylabel("Depth (m)")

    ###########################################################################
    # SMOOTHED UNCERTAINTY
    ###########################################################################
    
    pcm4 = axes[1, 1].pcolormesh(
        x_edges,
        z_edges,
        std_smooth,
        shading="auto",
        cmap=CMAP_STD,
        vmin=vmin_std,
        vmax=vmax_std,
    )
    
    axes[1, 1].set_title("Smoothed Vs Uncertainty")
    #axes[1, 1].set_xlabel("Distance along profile (m)")
    
    # Interface uncertainty

    ###########################################################################
    # COLOURBARS
    ###########################################################################
    
    cbar1 = fig.colorbar(
        pcm2,
        ax=axes[0, :],
        location="right"
    )
    cbar1.set_label("Vs (m/s)")
    
    cbar2 = fig.colorbar(
        pcm4,
        ax=axes[1, :],
        location="right"
    )
    cbar2.set_label("Vs Uncertainty (m/s)")

    ###########################################################################
    # FINISH
    ###########################################################################

    #plt.tight_layout()
    output_pdf = os.path.join(MODEL_DIR,f"{site}_{profile}_Vs_2D_Section.pdf")
    
    fig.savefig(
        output_pdf,
        format="pdf",
        bbox_inches="tight"
    )

    print(f"Saved: {output_pdf}")
    
    plt.show()
