import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import itertools

def moving_average(data, window_size):
    """
    Compute the moving average of a 1D array.

    Parameters:
        data (array-like): Input data series.
        window_size (int): Size of the moving window.

    Returns:
        np.ndarray: Array of smoothed values (same length as input).
    """
    data = np.asarray(data)
    if window_size <= 0 or window_size > len(data):
        raise ValueError("Invalid window size.")

    cumsum = np.cumsum(np.insert(data, 0, 0))  # pad with 0 for offset
    result = (cumsum[window_size:] - cumsum[:-window_size]) / window_size

    # Pad the beginning to return a full-length array (optional strategy)
    pad = np.full(window_size - 1, result[0])
    return np.concatenate((pad, result))

def calculate_msd(times, x_cm, lag_length):
    """
    Calculate mean square displacement (MSD) from a 1D trajectory with optional lag averaging.

    Parameters:
    - time: array-like
    - x_cm: array-like, position values over time.
    - max_lag: int or None, maximum lag (τ) to compute. Default is N/2.

    Returns:
    - lags: array of lag times (in same units as `time` if provided).
    - msd: numpy array of MSD values.
    """

    if isinstance(times, list):
        times = np.array(times)
    if isinstance(x_cm, list):
        x_cm = np.array(x_cm)

    N = x_cm.size
    dt = times[1] - times[0]
    msd = np.zeros(lag_length)
    lags = np.arange(1, lag_length + 1) * dt

    for tau in range(1, lag_length + 1):
        # Calculate squared displacements for each time lag
        # (x_cm[i + tau] - x_cm[i])**2 for all possible i
        displacements_squared = (x_cm[tau:] - x_cm[:-tau])**2
        msd[tau - 1] = np.mean(displacements_squared)

    return lags, msd

def wrap(cood, boxDim):
    rs = np.diff(list(itertools.combinations(cood, 2)), axis=1).squeeze()
    rs -= boxDim*np.floor(rs/boxDim+0.5) # PBC
    return rs

def animate_colormap(data, axs_labels=None, times=None, c_label=None, interval=50, sz=5, cm='bwr'):
    """
    Create an animation of a 2D colormap over time.

    This function takes a 3D numpy array representing time-dependent 2D data and 
    generates an animated colormap. The animation can be rendered in a Jupyter 
    notebook using the :class:`IPython.display.HTML` object.

    :param data: 
        A 3D array with shape (t, L, M), where t is the number of timesteps, 
        and L and M represent the dimensions of each data slice (e.g., rows and columns).
    :type data: numpy.ndarray
    :param axs_labels: 
        A list of length 2 containing strings for the x-axis and y-axis labels, 
        in the 0th and 1st positions respectively. Default is None.
    :type axs_labels: list of str, optional
    :param times: 
        A 1D array of length t, where each value represents the time corresponding 
        to each timestep. Default is None.
    :type times: numpy.ndarray, optional
    :param c_label: 
        A label for the colormap (color bar). Default is None.
    :type c_label: str, optional
    :param interval: 
        The delay in milliseconds between frames in the animation. Default is 50.
    :type interval: int, optional
    :param sz: 
        The size of the plot figure. Default is 5.
    :type sz: int, optional
    :param cm: 
        The colormap to use for the animation. Default is 'bwr'.
    :type cm: str, optional

    :return: 
        An animation object that can be rendered in a Jupyter notebook using 
        :class:`IPython.display.HTML`.
    :rtype: matplotlib.animation.FuncAnimation

    :examples:
        >>> import numpy as np
        >>> from IPython.display import HTML
        >>> data = np.random.random((50, 100, 100))  # Example 3D array (50 timesteps, 100x100 grid)
        >>> ani = animate_colormap(data, axs_labels=["X-axis", "Y-axis"], times=np.arange(50))
        >>> HTML(ani.to_jshtml())
    """
    def init():
        img.set_data(data[0])
        vmin = np.amin(data[0])
        vmax = np.amax(data[0])
        img.set_clim(vmin, vmax)
        if times is not None:
            ax.set(title=f"Time = {times[0]}")
        return (img,)

    def update(i):
        img.set_data(data[i])
        vmin = np.amin(data[i])
        vmax = np.amax(data[i])
        img.set_clim(vmin, vmax)
        if times is not None:
            ax.set(title=f"Time = {times[i]}")
        return (img,)

    fig, ax = plt.subplots(1, 1, figsize=(sz, sz))
    img = ax.imshow(data[0], cmap=cm, vmin=np.amin(data[0]), vmax=np.amax(data[0]))
    if axs_labels is not None:
        ax.set_xlabel(axs_labels[0])
        ax.set_ylabel(axs_labels[1])
    fig.colorbar(img, ax=ax, orientation="horizontal", label=c_label, pad=0.2)
    ani = animation.FuncAnimation(fig, update, frames=len(data), init_func=init, interval=interval, blit=True)
    plt.close()
    return ani