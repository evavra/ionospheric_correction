import os
import sys
import numpy as np
import xarray as xr
# import matplotlib.pyplot as plt\


def main():
    """
    Correct for discontinuities along seam boundaries and masks anomalous pixels.
    Shifts each subswath by the difference between median pixel values on either side of each seam.
    Left-most swath is used as reference.

    USAGE:
        python correct_subswath_local.py grd_file seam_file seam_width disp_max

    INPUT:
        grd_file   - path to grd file to correct
        seam_file  - textfile containing seam locations indices (left pixel, first pixel = 1)
        seam_width - number of pixels to average the difference between each subswath
        disp_max   - number of radians as a threshold to throw out

    OUTPUT:
        ph_correct - corrected grd file
    """

    if len(sys.argv) == 5:
        # Parse arguments
        grd_file   = sys.argv[1]
        seam_file  = sys.argv[2]
        seam_width = int(sys.argv[3])
        disp_max   = float(sys.argv[4])

        # Perform correction
        correct_subswath_local(grd_file, seam_file, seam_width, disp_max)
    else:
        # Return docstring if # arguments is wrong
        print(main.__doc__)
        sys.exit()

    return


def correct_subswath_local(grd_file, seam_file, seam_width, disp_max, file_name_corrected='ph_correct.grd'):
    # INPUT:
    # seam_width - number of pixels to average the difference between each subswath
    # disp_max   - number of radians as a threshold to throw out
    # noisy      - pixels whose amplitude > threshold (most tricky part)

    # boundary.txt saves the # of pixels for each subswath along range
    # can be generated when merge the swaths

    if os.path.exists(seam_file) is not True:
        print('The seam file does not exist!')
        sys.exit()
        return

    # Load data
    rng, azi, phase = read_grd(grd_file)
    seam_index      = np.loadtxt(seam_file, dtype=int)
    seam_index     -= 1 # Account for python indexing
    phase_corrected = np.copy(phase)

    L  = len(seam_index) # number of seams
    print(seam_index, phase.shape)
    subswaths = []


    for i in range(L):
        # Correct swath i+1
        # print(f'Correcting F{i + 2}...')

        # Get swath i + 1 columns
        if i == L - 1: # last swath
            start = seam_index[i]
            end   = phase.shape[1] 
        else:
            start = seam_index[i]
            end   = seam_index[i + 1]

        # Get median values columns around seam
        left_side  = phase[:, seam_index[i] - seam_width:seam_index[i]                 ].flatten()
        right_side = phase[:,              seam_index[i]:seam_index[i] + seam_width].flatten()
        offset     = np.nanmedian(left_side) - np.nanmedian(right_side)

        print()
        print('Seam:', seam_index[i] + 1)
        print('Left:', seam_index[i] - seam_width, seam_index[i])
        print('Right:', seam_index[i] + 1, seam_index[i] + seam_width + 1)
        phase_corrected[:, start:end] -= offset        

    # Set unwrapping errors to NaNs
    i_nans = np.abs(phase_corrected - np.nanmedian(phase_corrected)) > disp_max
    phase_corrected[i_nans] = np.nan

    # Write to grd file
    if os.path.exists(file_name_corrected):
        os.remove(file_name_corrected)
    
    write_grd(rng, azi, phase_corrected, file_name_corrected)

    return


def read_grd(file, flatten=False):
    """
    Read NetCDF grid file (x, y, z) to Numpy arrays.
    Specify flatten=True to get flattened arrays and original grid dimensions
        """
    
    with xr.open_dataset(file) as grd:
        try:
            x = grd.lon.values
            y = grd.lat.values
            z = grd.z.values

        except AttributeError:
            x = grd.x.values
            y = grd.y.values
            z = grd.z.values

    if flatten:
        dims = z.shape
        X, Y = np.meshgrid(x, y)
        X    = X.flatten()
        Y    = Y.flatten()
        Z    = z.flatten()

        return X, Y, Z, dims

    else: 
        return x, y, z


def write_grd(x, y, z, file_name, T=True, V=False):

    """
    Write gridded data to GMT-compatible NetCDF file. 
    Requires GMT installation to modify grid registration (-T option).
    """

    data = xr.Dataset({'z': (['y', 'x'], z)}, coords={'y':(['y'], y), 'x':(['x'], x)})

    # Fix headers
    x_inc = np.min(np.diff(data.x.values))
    y_inc = np.min(np.diff(data.y.values))

    x_min_orig = np.min(data.x.values)
    x_max_orig = np.max(data.x.values)
    y_min_orig = np.min(data.y.values)
    y_max_orig = np.max(data.y.values)

    data.x.attrs['actual_range'] = np.array([x_min_orig, x_max_orig])
    data.y.attrs['actual_range'] = np.array([y_min_orig, y_max_orig])

    data.to_netcdf(path=file_name, mode='w')

    if T:
        x_min = np.min(data.x.values) - x_inc/2 
        x_max = np.max(data.x.values) + x_inc/2
        y_min = np.min(data.y.values) - y_inc/2
        y_max = np.max(data.y.values) + y_inc/2

        bounds = ' -R' + str(x_min) + '/' + str(x_max) + '/' + str(y_min) + '/' + str(y_max)
        cmd = 'gmt grdedit ' + file_name + ' ' + bounds + ' -T' 
        if V:
            print(cmd)
        os.system(cmd)

    return


if __name__ == '__main__':
    main()


