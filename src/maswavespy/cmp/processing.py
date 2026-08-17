def process_cmp(cmp_loc, config):
    from sys import platform as sys_pf
    if sys_pf == 'darwin':
        import matplotlib
        matplotlib.use("TkAgg")
    from maswavespy import wavefield
    from maswavespy import inversion
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import pickle
    import os
    print(os.getcwd())
    # Input filenames
    cmp_dir = (config["cmp_sorting"]["output_dir"])
    file_name = (f"{cmp_dir}/CMP_{cmp_loc:.2f}.dat")
    geometry_file = (f"{cmp_dir}/CMP_{cmp_loc:.2f}_geometry.csv")
    # site and profile definitions
    site = config["project"]["site"]
    profile = config["project"]["profile"]
    #Sampling Rates
    fs = config["processing"]["sampling_rate"]
    f_pick_min = (config["processing"]["f_pick_min"])
    # Dispersion calculation settings
    disp = config["processing"]["dispersion"]
    cT_min = disp["cT_min"]
    cT_max = disp["cT_max"]
    cT_step = disp["cT_step"]
    # INFO from yaml file
    inventory = pd.read_csv(
        os.path.join(
            config["cmp_sorting"]["output_dir"],
            "CMP_inventory.csv"
        )
    )

    row = inventory.loc[
        np.isclose(inventory["x1"], cmp_loc)
    ].iloc[0]
    n = int(row["n"])
    header_lines = int(row["header_lines"])
    direction = row["direction"]
    dx = float(row["dx"])
    x1 = float(row["x1"])
    
    # Importing the data into rec_cmp - data from textfile as written out as CMP
    rec_cmp = wavefield.RecordMC.import_from_textfile(
        site=site,
        profile=profile,
        file_name=file_name,
        header_lines=header_lines,
        n=n,
        direction=direction,
        dx=dx,
        x1=x1,
        fs=fs,
        f_pick_min=f_pick_min,
        geometry_file=geometry_file
    )
    """
    rec_cmp = wavefield.RecordMC.import_from_textfile(
        site=site,
        profile=profile,
        file_name=file_name,
        header_lines=5,
        direction='cmp',
        dx=dx,
        x1=x1,
        fs=fs,
        f_pick_min=f_pick_min,
        geometry_file=geometry_file
    )
    """
    # Some output to check data consistency
    print(f"Offsets:{rec_cmp.offsets}") # offsets in CMP
    print(f"Sources: {rec_cmp.sources}") # Sources are actually CMP locations
    print(f"Receivers: {rec_cmp.receivers}") # Geophone locations 
    
    # Some more output
    print(f"Direction: {rec_cmp.direction}")
    print(f"Offsets: {rec_cmp.offsets[:10]}")
    
    # Plotting the data
    rec_cmp.plot_data(du=0.75,normalized=True,filled=True)
    plt.close()
    # Calculating Dispersion
    print("Calculating Dispersion")
    
    # creating a dispersion curve element
    edc_cmp = rec_cmp.element_dc(cT_min,cT_max,cT_step)

    # plot the dispersion image
    #edc_cmp.plot_dispersion_image(f_min,f_max)

    # Pick Dispersion
    # Launch a GUI for interactive pickint
    disp = config["processing"]["dispersion"]
    f_min = config["processing"]["dispersion"]["f_min"]
    f_max = disp["f_max"]
    model_dir = config["output"]["model_dir"]

    os.makedirs(
        model_dir,
        exist_ok=True
    )
    print(f"Launching GUI for dispersion curve picking for CMP {cmp_loc:.2f}")
    edc_cmp.pick_dc(f_min,f_max)

    # get the dispersion curve infor
    dc_dict = edc_cmp.get_inversion_dict()

    pick_file = (f"{model_dir}/"f"CMP_{cmp_loc:.2f}_pick.pkl")

    with open(pick_file, "wb") as f:
        pickle.dump(
            dc_dict,
            f
        )
    # Import initial soil model parameters from csv file
    filename_initial = (config["input"]["initial_model"])
    initial_parameters = pd.read_csv(filename_initial)
    
    h = np.array(initial_parameters['h [m]'].values[0:-1], dtype='float64')
    n = int(len(h))
    Vs = np.array(initial_parameters['Vs [m/s]'].values, dtype='float64')
    rho = np.array(initial_parameters['rho [kg/m3]'].values, dtype='float64')
    Vp = []
    n_unsat = 0; nu = None
    for item in range(len(initial_parameters['saturated/unsaturated'].values)):
        if initial_parameters['saturated/unsaturated'].values[item] == 'unsat':
            nu = initial_parameters['nu [-]'].values[item]
            Vp.append(np.sqrt((2*(1-nu))/(1-2*nu))*Vs[item])
            n_unsat = n_unsat + 1
        else:
            Vp.append(initial_parameters['Vp [m/s]'].values[item])
    Vp = np.array(Vp, dtype='float64')
    # Print message to user
    print('The sample dispersion curve has been imported.')
    print('The initial soil model parameters have been imported.')

    # Getting the parameters from the dispersion dictionary
    c_mean = dc_dict["c_mean"]
    c_low = dc_dict["c_low"]
    c_up = dc_dict["c_up"]
    wavelengths = dc_dict["wavelength"]

    # Initialize an inversion object.
    
    inv_TestSite = inversion.InvertDC(site, profile, c_mean, c_low, c_up, wavelengths)
    # Print message to user
    print('An inversion (InvertDC) object has been initialized.') 
    
    # Initialize the inversion routine. The inversion is conducted using 
    # Monte Carlo sampling as described in Olafsdottir et al. (2020).
    
    # Range for testing phase velocity
    ct = (config["processing"]["inversion"]["c_test"])

    c_test = {"min": ct["min"],"max": ct["max"],"step": ct["step"],"delta_c": ct["delta_c"]}

    # Initial model parameters
    initial = {'n' : n,
               'n_unsat' : n_unsat, # that was n_unsat before
               'alpha' : Vp,
               'nu_unsat' : 0.35,
               'alpha_sat' : 1500,
               'beta' : Vs,
               'rho' : rho,
               'h' : h,
               'reversals' : 0}      

    # Inversion algorithm settings. See further in Olafsdottir et al. (2020).
    setcfg = (config["processing"]["inversion"]["settings"])

    settings = {"run": setcfg["run"],
                "bs": setcfg["bs"],
                "bh": setcfg["bh"],
                "N_max": setcfg["N_max"],
                }
    max_depth = (config["processing"]["max_depth"])
    # View the initial shear wave velocity profile.
    # Compute the associated dispersion curve and show relative to the
    # experimental data. The misfit value is printed to the screen.
    inv_TestSite.view_initial(initial, max_depth, c_test, col='crimson',
                              DC_yaxis='linear', fig=None, ax=None,
                              figwidth=16, figheight=12, return_ct=False)

    # Print message to user
    print('The initial estimate of the Vs profile and the corresponding theoretical DC have been plotted.')
    plt.show()
    # Start the inversion analysis (optimization) process.
    print('Inversion initiated.')
    inv_TestSite.mc_inversion(c_test, initial, settings)
    
    # Plot sampled Vs profiles and associated dispersion curves
    inv_TestSite.plot_sampled(max_depth, runs='all', figwidth=16,
                              figheight=12, col_map='viridis', 
                              colorbar=True, DC_yaxis='linear',
                              return_axes=False, show_exp_dc=True)

    # Print message to user
    print('All runs completed.')
    print('The sampled Vs profiles and the corresponding theoretical DCs have been plotted.')
    plt.show()
    # Pickle the inversion object
    file = f'{site}_{profile}_CMP{cmp_loc}.00_inversion'
    inv_TestSite.save_to_pickle(file)
    # Print message to user
    print('The InvertDC object has been saved to disk as ' + file + '.p using pickle.')
    
    # Plot sampled Vs profiles whose associated dispersion curves fall within the
    # boundaries defined by c_low and c_up at all wavelengths
    inv_TestSite.plot_within_boundaries(max_depth, show_all=True,
                                        runs='all', figwidth=16, figheight=12, 
                                        col_map='viridis', colorbar=True,
                                        DC_yaxis='linear', return_axes=False)

    # Print message to user
    print('The set of accepted Vs profiles and the corresponding theoretical DCs have been plotted.')
    plt.show()

    # #### Post-processing 
    # 
    # #### Compute and plot the median shear wave velocity profile (defined
    # in terms of shear wave velocity and depth of layer interfaces) and
    # the 90-th percentiles of each parameter. The associated theoretical
    # dispersion curve is also computed and shown relative to the experimetnal
    # data. The mean shear wave velocity profile can be obtained in a
    # comparable way using the inv_TestSite.mean_profile method. (See further
    # in the documentation of inversion.py.)
    percentiles = [10,90]
    TestSite_median_profile = inv_TestSite.median_profile(q=percentiles,
                                                          dataset='selected')
    len(TestSite_median_profile)
    #print(TestSite_median_profile)
    fig, ax = inv_TestSite.plot_profile(TestSite_median_profile, max_depth,
                                        c_test, initial, col='red',
                                        up_low=True, fig=None, ax=None, 
                                        return_axes=True, return_ct=False)
    # Print message to user
    print('The median of accepted Vs profiles has been computed.')
    print('The median profile and the corresponding theoretical DC have been plotted.')
    plt.show()
    # Now plot some mean profiles
    TestSite_mean_profile = inv_TestSite.mean_profile(stdev=True, no_stdev=1, dataset='selected')
    """
    print(TestSite_mean_profile['beta'])
    fig, ax = inv_TestSite.plot_profile(TestSite_mean_profile, max_depth,
                                        c_test, initial, col='red',
                                        up_low=True, fig=None, ax=None,
                                        return_axes=True, return_ct=False)
    print(TestSite_mean_profile)
    """
    # writing out the file
    model_dir = (config["output"]["model_dir"])

    os.makedirs(model_dir,exist_ok=True)

    file = (f"{model_dir}/"f"__CMP_{cmp_loc:.2f}_mean")

    with open(f'{file}.pkl', 'wb') as f:
        pickle.dump(TestSite_mean_profile, f)
        # Print message to user
        print('The mean velocity and standard deviations have been saved to disk as ' + file + '.p using pickle.')

        # Look at the lower misfit profiles
        lowest_misfit_profiles = {}
        no_profiles = 10
        # Ensure that at least no_profiles fall within the experimental DC boundaries
        no_profiles_checked = min(no_profiles, len(inv_TestSite.selected['beta']))
        for no in range(-1*no_profiles_checked,0):
            profile_dict = {'beta': inv_TestSite.selected['beta'][no],
                            'z': inv_TestSite.selected['z'][no]}
            if no == -1*no_profiles_checked:
                fig, ax = inv_TestSite.plot_profile(profile_dict, max_depth,
                                                    c_test, initial,
                                                    col='gray', up_low=False,
                                                    DC_yaxis='linear',
                                                    fig=None, ax=None,
                                                    return_axes=True,
                                                    show_legend=True)
            else:
                inv_TestSite.plot_profile(profile_dict, max_depth, c_test,
                                          initial, col='gray', up_low=False,
                                          DC_yaxis='linear', fig=fig, ax=ax,
                                          show_legend=False)
                lowest_misfit_profiles[no] = profile_dict

        # Print message to user
        print('The ' + str(no_profiles) +
              ' lowest-misfit Vs profiles have been identified.')
        print('The set of lowest-misfit profiles and the corresponding theoretical DCs have been plotted')


        # #### Post-processing 
        # 
        # #### Compute the average shear wave velocity (Vsz) for the top most z=5 m, z=10 m, z=20 m and z=30 m using (i) the median Vs profile and (ii) the lowest-misfit Vs profile.
        # 
        depth = [5, 10, 20, 30]
        layer_parameter = 'z'

        # Median Vs profile
        Vsz_median = inv_TestSite.compute_vsz(depth,
                                              TestSite_median_profile['beta'], 
                                              TestSite_median_profile['z'],
                                              layer_parameter)
        print('Median Vs profile, z and Vsz values')
        print(Vsz_median[0]) # Depths (z)
        print([round(val, 2) for val in Vsz_median[1]]) # Computed Vsz values

        # Lowest-misfit Vs profile
        Vsz_lowest_misfit = inv_TestSite.compute_vsz(depth,inv_TestSite.selected['beta'][-1], inv_TestSite.selected['z'][-1], layer_parameter)
        print('Lowest-misfit Vs profile, z and Vsz values')
        print(Vsz_lowest_misfit[0]) # Depths (z)
        print([round(val, 2) for val in Vsz_lowest_misfit[1]]) # Computed Vsz values
