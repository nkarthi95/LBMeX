import yt
import numpy as np
import csv

def read_amrex_data(ds, boxDim):
    """
    Reads and organizes data from AMReX plotfiles.

    This function parses and reshapes data created by AMReX when the output files are in the native plotfile format. The plotfiles are pre-processed using `yt` to create a data series object (`ds`), which is then used to retrieve and reshape the data into the desired format.

    Parameters
    ----------
    ds : yt.data_objects.static_output.Dataset
        The data series object created by `yt` from AMReX plotfiles.
    boxDim : list or numpy.ndarray
        The dimensions of the simulation box, provided as a list or NumPy array of size 2 (for 2D) or 3 (for 3D).
    begin : int, optional
        The starting index of the data fields to retrieve. Defaults to 0.
    end : int, optional
        The ending index (exclusive) of the data fields to retrieve. Defaults to 5.

    Returns
    -------
    numpy.ndarray
        A NumPy array containing the organized data. The shape of the array is:
        - `(M, nx, ny, nz)` for 3D data, where `M = end - begin` and `nx`, `ny`, `nz` are the dimensions of the box.
        - `(M, nx, ny)` for 2D data, where `M = end - begin` and `nx`, `ny` are the dimensions of the box.

    Notes
    -----
    - The input `boxDim` is converted to a NumPy array internally if provided as a list.
    - Data fields from the `ds` object are accessed sequentially, reshaped according to `boxDim`, and stored in the output array.
    - The function uses `yt`'s `to_dataframe` method for sorting and reshaping the data.

    Examples
    --------
    Process a dataset from AMReX plotfiles and reshape it for analysis:

    >>> import numpy as np
    >>> import yt
    >>> ds = yt.load("plt00000")  # Pre-processed AMReX plotfile
    >>> boxDim = [64, 64, 64]  # Dimensions of the simulation box
    >>> data = read_amrex_data(ds, boxDim, begin=0, end=3)
    >>> print(data.shape)
    (3, 64, 64, 64)
    """
    if type(boxDim) == list:
        boxDim = np.array(boxDim)

    data_fields = ds.field_list
    all_data_level_0 = ds.covering_grid(level=0, left_edge=[0, 0.0, 0.0], dims=boxDim)

    run_data = {}

    for i, field_info in enumerate(data_fields):
        field = field_info[-1]
        run_data[field] = np.array(all_data_level_0[field])
    
    return run_data

def read_csv(path):
    ls = []
    with open(path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=" ", quotechar="|")
        
        for idx, row in enumerate(reader):
            if not row or not row[0].strip():
                continue  # Skip empty rows

            tokens = row[0].strip().split(",")[:-1]  # Remove trailing empty after final comma

            # Attempt to convert all to float
            try:
                parsed = [float(tok) for tok in tokens]
            except ValueError:
                if idx == 0:
                    # First row failed to convert — assume it's a header
                    continue
                else:
                    # Data row failed to convert — preserve as strings
                    parsed = tokens

            ls.append(parsed)

    return np.array(ls, dtype=object if any(isinstance(i, str) for row in ls for i in row) else float)