import os
from collections import defaultdict

import numpy as np
import pandas as pd
from obspy import read

def create_cmp_gathers(
        input_csv,
        output_dir,
        cmp_bin_factor=2,
        forward_only=True):

    """
    Generate CMP gathers from SEG2 data.

    Parameters
    ----------
    cmp_bin_factor : float
        CMP bin size multiplier relative to receiver spacing.
        Example: 2 -> CMP bins = 2 * receiver_spacing.

    input_csv : str
        CSV file containing SEG2 filenames in the first column.

    output_dir : str
        Directory where CMP gathers and geometry files are written.

    forward_only : bool, optional
        Keep only forward shots (source < receiver).
        Default is True.

    Returns
    -------
    inventory_df : pandas.DataFrame
        Inventory table written to CMP_inventory.csv.

    cmp_gathers : dict
        Dictionary containing all CMP gathers.
    """

    os.makedirs(output_dir, exist_ok=True)

    # ------------------------------------------------------
    # Read file list
    # ------------------------------------------------------

    df = pd.read_csv(input_csv)
    seg2_files = df.iloc[:, 0].tolist()

    # ------------------------------------------------------
    # Storage
    # ------------------------------------------------------

    cmp_gathers = defaultdict(list)
    inventory = []

    acq_date = None
    acq_time = None
    sampling_rate = None
    receiver_spacing = None

    # ------------------------------------------------------
    # Read SEG2 files
    # ------------------------------------------------------

    for segfile in seg2_files:

        print(f"Reading {segfile}")

        st = read(segfile)

        # --------------------------------------------------
        # Common acquisition information
        # --------------------------------------------------

        if acq_date is None:

            hdr = st[0].stats.seg2

            acq_date = hdr["ACQUISITION_DATE"]
            acq_time = hdr["ACQUISITION_TIME"]

            delay = float(hdr["DELAY"])
            sampling_rate = 1.0 / st[0].stats.delta

            r1 = float(st[0].stats.seg2["RECEIVER_LOCATION"])
            r2 = float(st[1].stats.seg2["RECEIVER_LOCATION"])

            receiver_spacing = abs(r2 - r1)

            cmp_bin_size = receiver_spacing * cmp_bin_factor

            print(f"Sampling rate: {sampling_rate:.1f} Hz")
            print(f"Receiver spacing: {receiver_spacing:.2f} m")
            print(f"Shot delay: {delay:.4f} s")
            print(f"CMP bin size: {cmp_bin_size:.2f} m")

        # --------------------------------------------------
        # Loop through traces
        # --------------------------------------------------

        for tr in st:

            hdr = tr.stats.seg2

            source_x = float(hdr["SOURCE_LOCATION"])
            receiver_x = float(hdr["RECEIVER_LOCATION"])

            cmp_x = 0.5 * (source_x + receiver_x)

            cmp_bin = (
                round(cmp_x / cmp_bin_size)
                * cmp_bin_size
            )

            offset = receiver_x - source_x

            # ----------------------------------------------
            # Forward-only option
            # ----------------------------------------------

            if forward_only:

                if source_x >= receiver_x:
                    continue

                geometry_source = cmp_bin
                geometry_receiver = offset

            else:

                geometry_source = source_x
                geometry_receiver = receiver_x

            dt = tr.stats.delta
            delay = float(hdr["DELAY"])

            trace_data = tr.data.astype(np.float64)

            if delay < 0:

                nshift = int(round(abs(delay) / dt))

                if nshift < len(trace_data):
                    trace_data = trace_data[nshift:]
                else:
                    continue

            elif delay > 0:

                nshift = int(round(delay / dt))

                trace_data = np.concatenate(
                    [
                        np.zeros(
                            nshift,
                            dtype=trace_data.dtype,
                        ),
                        trace_data,
                    ]
                )

            cmp_gathers[cmp_bin].append(
                {
                    "trace": trace_data,
                    "source": geometry_source,
                    "receiver": geometry_receiver,
                    "offset": offset,
                    "shot": hdr["SHOT_SEQUENCE_NUMBER"],
                    "file": segfile,
                }
            )

    # ------------------------------------------------------
    # Write CMP gathers
    # ------------------------------------------------------

    print("\nWriting CMP gathers\n")

    record_counter = 1

    for cmp_loc in sorted(cmp_gathers.keys()):

        traces = cmp_gathers[cmp_loc]

        traces.sort(key=lambda x: x["offset"])

        min_length = min(
            len(t["trace"])
            for t in traces
        )

        data = np.column_stack(
            [
                t["trace"][:min_length]
                for t in traces
            ]
        )

        ntr = data.shape[1]

        waveform_file = os.path.join(
            output_dir,
            f"CMP_{cmp_loc:.2f}.dat",
        )

        with open(waveform_file, "w") as f:

            f.write(f"{cmp_loc:.2f}\n")

            f.write(
                f"{acq_date}\t"
                f"{acq_time}\t"
                f"{sampling_rate:.0f}\n"
            )

            f.write(
                f"{receiver_spacing:.2f}\t"
                f"{cmp_loc:.2f}\n"
            )

            f.write("\n")

            channel_line = "\t".join(
                [
                    f"Channel {i + 1}"
                    for i in range(ntr)
                ]
            )

            f.write(channel_line + "\n")

            np.savetxt(
                f,
                data,
                delimiter="\t",
                fmt="%.15g",
            )

        geometry_file = os.path.join(
            output_dir,
            f"CMP_{cmp_loc:.2f}_geometry.csv",
        )

        geometry = pd.DataFrame(
            {
                "trace": np.arange(1, ntr + 1),
                "source": [t["source"] for t in traces],
                "receiver": [t["receiver"] for t in traces],
                "offset": [t["offset"] for t in traces],
                "cmp": np.full(ntr, cmp_loc),
                "file": [
                    os.path.basename(t["file"])
                    for t in traces
                ],
                "shot": [t["shot"] for t in traces],
            }
        )

        geometry.to_csv(
            geometry_file,
            index=False,
        )

        inventory.append(
            {
                "record_id": f"r{record_counter}",
                "site": "CMP",
                "profile": "CMP",
                "file_name": os.path.basename(
                    waveform_file
                ),
                "header_lines": 5,
                "n": ntr,
                "direction": "cmp",
                "dx": receiver_spacing,
                "x1": cmp_loc,
                "geometry_file": os.path.basename(
                    geometry_file
                ),
            }
        )

        print(
            f"CMP {cmp_loc:.2f} m : "
            f"{ntr} traces"
        )

        record_counter += 1

    # ------------------------------------------------------
    # Write inventory
    # ------------------------------------------------------

    inventory_df = pd.DataFrame(inventory)

    inventory_file = os.path.join(
        output_dir,
        "CMP_inventory.csv",
    )

    inventory_df.to_csv(
        inventory_file,
        index=False,
    )

    print("\nFinished.")
    print(f"Inventory written to: {inventory_file}")
    print(f"CMP gathers written to: {output_dir}")

    inventory_file = (
        f"{output_dir}/CMP_inventory.csv"
    )

    inventory_df.to_csv(
        inventory_file,
        index=False
    )
    return sorted(cmp_gathers.keys())
    """
    print(inventory_df.columns)
    print(inventory_df.head())
    cmp_locations = sorted(
        inventory_df["x1"]
        .unique()
        .tolist()
    )

    return cmp_locations
    """
