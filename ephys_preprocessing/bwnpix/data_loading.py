import numpy as np
import sys
import os
import glob

# Append the path to the brainwide-npix/bwnpix directory to the system path
sys.path.append(os.path.join(os.path.dirname(__file__), '../brainwide-npix/bwnpix'))
from load import *
from multiprobe import *
from bwnpix.multiprobe import MultiprobeEphysExperiment
from scipy.stats import zscore
import scipy
import tqdm
from pylab import *

import pandas as pd
from matplotlib.colors import ListedColormap
import scipy.io
# Define paths for data storage and retrieval
cache_root = '/path/to/cache_root/'
behavior_path = '/path/to/npix_behavior/'
fig_path = '/path/to/figures/'
data_export_path = '/path/to/export/'
ephys_path = '/path/to/ephys/localProcessedData/'

def get_good_recordings(condition):
    '''
    Returns a list of recordings for a specific experimental condition.

    Args:
        condition (str): The experimental condition. Should be one of: 
                         'both', 'hungry', 'thirsty', 'sated'.

    Returns:
        list: A list of recording names that match the specified condition.
    '''
    jackson = [ 
            "W01_20250716",
            #INSERT NAME OF RECORDING HERE
        ]

    if condition=="both":
        return both
    elif condition=="hungry":
        return hungry
    elif condition=="thirsty":
        return thirsty
    elif condition=="sated":
        return sated
    elif condition == "jackson":
        return jackson
    
def cache_data(recordings):
    '''
    Caches MultiprobeData objects for a list of recordings.

    Args:
        recordings (list): A list of recording names in the format 'subject_date', 
                           e.g., ['mon3_20231212', 'mon3_20231211'].

    Returns:
        list: A list of paths to the cached .pkl files.
    '''
    from bwnpix import load
    from bwnpix import anatomy

    datasets = []
    for recording in recordings:
        mname, mdate = recording.split('_')
        path = os.path.join(ephys_path, f"{mname}/catgt_{recording}_g0")
        print('loading' + path)
        if not os.path.exists(path):
            raise ValueError
        else:
            datasets.append([mname, mdate, path])
    
    use_cache = False # Set to true after this has been run once
    tract_names = [""] # tract name for easy reference
    all_paths = []

    for dataset in tqdm.tqdm(datasets):
        pathname = os.path.join(data_export_path, f"{'_'.join(dataset[:2])}.pkl")
        all_paths.append(pathname)
        if (not use_cache) or (not os.path.exists(pathname)):
            # if using bombcell for unit quality, set qc='bombcell', otherwise 'ks'
            mpd = load.MultiprobeData(dataset[0], dataset[1], qc='bombcell', load_lfp=False)
            for j in range(1): # add 4 probes
                mpd.add_experiment(os.path.join(dataset[2], f"{dataset[0]}_{dataset[1]}_g0_imec{j}",f"imec{j}_ks4"), tract_names[j])
            mpd.combine_experiments()
            try:
                event_map = {"xa_1": "trigger"}
                mpd.add_events(dataset[2], event_map)
            except:
                try:
                    event_map = {"xa_2": "trigger"}
                    mpd.add_events(dataset[2], event_map)
                except:
                    pass
            mpd.save_data(pathname)
    return all_paths

# Set numpy print options
np.set_printoptions(threshold=100)

## COLORS 
# Sets the dictionaries for the colors of behavior (reorient vs approach) and reward (food vs water)
behavior_colormap = ListedColormap(["#009AB1", "#9900DF", "#005F1E"])
reward_colormap = ListedColormap(["#CC470A", "#1A63A6"])

behavior_map = {'Food approach': 0,
                'Water approach': 0,
                'Food reorient left': 1,
                'Water reorient left': 1,
                'Food reorient right': 2,
                'Water reorient right': 2}

reward_map = {'Food approach': 0,
                'Water approach': 1,
                'Food reorient left': 0,
                'Water reorient left': 1,
                'Food reorient right': 0,
                'Water reorient right': 1}

def multiple_tests(pvals):
    '''
    Performs multiple hypothesis testing correction using the Benjamini/Hochberg method.

    Args:
        pvals (np.ndarray): A 1D numpy array of p-values that need correction.

    Returns:
        np.ndarray: A 1D numpy array of corrected p-values after applying the 
                    false discovery rate (FDR) control.
    '''
    from statsmodels.stats.multitest import multipletests  # Importing multipletests from statsmodels for p-value correction
    # Perform the multiple tests correction using Benjamini/Hochberg method (FDR)
    _, new_pvals, _, _ = multipletests(np.nan_to_num(pvals), method="fdr_bh")  
    return new_pvals  # Return the corrected p-values

def smooth_and_cache_fr(all_paths, smooth_size=5, save_dir=data_export_path):
    '''
    Smooths and caches firing rate data for multiple recordings.

    Args:
        all_paths (list): A list of paths to .pkl files containing MultiprobeEphysExperiment data.
        smooth_size (int, optional): The size of the smoothing window. Defaults to 5.
        save_dir (str, optional): The directory to save the smoothed firing rate data. 
                                  Defaults to the global variable `data_export_path`.
    '''
    def smooth_data(data, bin_t=0.01, smooth=5):
        fr, bins = data.firing_rates_over_recording(bin_size=bin_t, smooth_size=smooth)
        zfr = zscore(fr, 1)
        bins_fr = bins
        return fr, zfr 

    for path in tqdm.tqdm(all_paths):
        data = MultiprobeEphysExperiment()
        data.load_units_from_ks(path, load_lfp=False)
        # Get filename from path
        dset = os.path.basename(path).split('.')[0].split('_')
        data._mouse_name = dset[0]; data._recording = f"{dset[0]}_{dset[1]}"
        res = smooth_data(data, smooth=smooth_size)
        fr = res[0]
        zfr = res[1]
        np.save(os.path.join(save_dir, f"{dset[0]}_{dset[1]}_fr_{smooth_size}.npy"), fr)
        np.save(os.path.join(save_dir, f"{dset[0]}_{dset[1]}_zfr_{smooth_size}.npy"), zfr)

def load_behavior(recording_name):
    '''
    Loads behavior data for a given recording.

    Args:
        recording_name (str): The name of the recording.

    Returns:
        dict: A dictionary containing the loaded behavior data.
    '''
    behavior = scipy.io.loadmat(behavior_path + recording_name + '.mat')
    return behavior

def convert_cam_time(frame_number, sync_frames, npix_sync_times):
    '''
    Converts a frame number from the camera's clock to Neuropixels time.

    Args:
        frame_number (int): The frame number from the camera.
        sync_frames (array-like): An array of frame numbers corresponding to synchronization events.
        npix_sync_times (array-like): An array of Neuropixels timestamps corresponding to synchronization events.

    Returns:
        float: The Neuropixels time corresponding to the input frame number.
    '''
    # Find the sync frame that is closest to the frame number
    sync_frame_idx = np.max([0, np.searchsorted(sync_frames, frame_number, side='right')-1])
    try:
        npix_sync_time = npix_sync_times[sync_frame_idx]
    except:
        sync_frame_idx = len(npix_sync_times)-1
        npix_sync_time = npix_sync_times[sync_frame_idx]
    # Multiple cases to avoid issues at start or end of recording
    try:
        time_between_syncs = npix_sync_times[sync_frame_idx+1] - npix_sync_times[sync_frame_idx]
        framerate = time_between_syncs / (sync_frames[sync_frame_idx+1] - sync_frames[sync_frame_idx])
    except:
        time_between_syncs = npix_sync_times[sync_frame_idx] - npix_sync_times[sync_frame_idx-1]
        framerate = time_between_syncs / (sync_frames[sync_frame_idx] - sync_frames[sync_frame_idx-1])
    # Find number of elapsed frames
    frames_elapsed = frame_number - sync_frames[sync_frame_idx]
    # Add elapsed time to sync time, on npix clock
    npix_time = npix_sync_time + (frames_elapsed * framerate)
    return npix_time

def convert_computer_time(computer_time, computer_sync_times, npix_sync_times):
    '''
    Converts computer time to Neuropixels time in seconds.

    Args:
        computer_time (float): The time in seconds from the computer's clock.
        computer_sync_times (array-like): An array of computer timestamps corresponding to synchronization events.
        npix_sync_times (array-like): An array of Neuropixels timestamps corresponding to synchronization events.

    Returns:
        float: The Neuropixels time corresponding to the input computer time.
    '''
    # Find computer sync that is closest to input time
    sync_idx = np.max([0, np.searchsorted(computer_sync_times, computer_time, side='right')-1])
    npix_sync_time = npix_sync_times[sync_idx]
    time_elapsed = computer_time - computer_sync_times[sync_idx]
    npix_time = npix_sync_time + time_elapsed
    return npix_time

def get_data_by_windows(windows, data, bin_t=0.01):
    '''Get firing rates for bout windows'''
    d = []
    for start, end in windows:
        d.append(data[:,int(start/bin_t):int(end/bin_t)])
    return d

def compute_cd_cov(x,y):
    """
    Inputs:
        x: Ncell x Ntrial 
    """
    v = (x.mean(1)-y.mean(1))/np.sqrt(np.var(x) + np.var(y))
    return v/np.sum(np.abs(v))

def load_recording_data(recording_names, atlas, dir=data_export_path):
    '''
    Loads recording data for a list of recording names.

    Args:
        recording_names (list): A list of recording names.
        atlas (object): An atlas object for histology data.
        dir (str, optional): The directory where the data is stored. 
                             Defaults to 'Z:/motiv/data/export/'.

    Returns:
        list: A list of MultiprobeEphysExperiment objects containing the loaded data.
    '''
    all_data = []
    for recording_name in tqdm.tqdm(recording_names):
        path = data_export_path+recording_name+'.pkl'
        print(path)
        print(f"Loading recording: {recording_name}")
        data = MultiprobeEphysExperiment()
        data.load_units_from_ks(path, load_lfp=False)
        data._mouse_name = recording_name.split('_')[0]
        data._recording = recording_name
        data.load_histology(atlas)
        all_data.append(data)
    return all_data

def get_fr_for_recordings(recording_names, zscored=True, smooth=5, dir=data_export_path):
    '''
    Retrieves firing rate data for a list of recordings.

    Args:
        recording_names (list): A list of recording names.
        zscored (bool, optional): Whether to retrieve z-scored firing rates. Defaults to True.
        smooth (int, optional): The smoothing window size used for the firing rate calculation. Defaults to 5.
        dir (str, optional): The directory where the firing rate data is stored. 
                             Defaults to 'Z:/motiv/data/export/'.

    Returns:
        list: A list of NumPy arrays containing the firing rate data for each recording.
    '''
    all_fr = []  # Initialize an empty list to hold firing rate data for all recordings
    for recording_name in recording_names:  # Iterate over each recording name
        # Construct the file path based on whether z-score is requested
        if zscored:
            path = glob.glob(os.path.join(dir, f'{recording_name}_zfr_{smooth}.npy'))[0]  # Find z-scored firing rate file
        else:
            path = glob.glob(os.path.join(dir, f'{recording_name}_fr_{smooth}.npy'))[0]  # Find normal firing rate file

        print(f"Loading recording fr: {recording_name}")  # Print the name of the recording being loaded
        all_fr.append(np.load(path, mmap_mode='r'))  # Load the firing rate data and append it to the list
    return all_fr  # Return the list of firing rate data

def get_annotated_syllables(recording_name, daq_sync=None, syllable_type="supervised"): 
    '''Get behavior syllable data for a recording.
    Daq_sync is the nidaq-recorded event times for each sync, on neuropixels time.
    Syllables are assigned a unique index.
    Return a pandas dataframe for each syllable bout containing:
    - syllable bout index
    - start frame
    - start time (neuropixels time)
    - end frame
    - end time (neuropixels time)
    - duration (in seconds)
    - syllable name
    - preceeding reward type
    - time from previous reward
    - subsequent reward type
    - time to next reward
    - preceding syllable name
    - subsequent syllable name
    - position at start of syllable
    - position at end of syllable
    '''
    mat = load_behavior(recording_name)
    bouts = []
    sync_frames_foraging = mat['locs']['locfz'][0][0].flatten()-1
    sync_times_computer = mat['syncs'].flatten()

    if syllable_type == "supervised":
        syllable_name_to_index = {
            "None": 0,
            "consume": 1,
            "reorient left": 2,
            "approach": 3, # approach
            "reorient right": 4,
            "turn": 5  
        }
        syllable_by_frame = mat['behave']['behavior'][0][0].flatten()

    else:
        syllable_name_to_index = {
        "approach": 2,
        "consume": 1,
        "reorient left": 3,
        "reorient right": 5,
        "curled pause/grooming": 4,
        "investigate corner": 6,
        "pause/turn & scratch": 7,
        "head nod": 8,
        "rear/investigate/nod": 9
        }
        syllable_by_frame = mat['kpmsfz']['syllreindexed'][0][0].flatten()

    syllable_index_to_name = {v: k for k, v in syllable_name_to_index.items()}

    foods = mat['rewcol']['food'][0][0].flatten()
    waters = mat['rewcol']['water'][0][0].flatten()
    outcome = np.array(len(foods)*['food'] + len(waters)*['water'])
    outcomes = outcome[np.argsort(np.concatenate([foods, waters]))]
    total_water = np.cumsum(1.0*(outcomes=='water'))
    total_food = np.cumsum(1.0*(outcomes=='food'))
    reward_times = np.sort(np.concatenate([foods, waters]))
    # convert reward times to npix time
    reward_times = np.array([convert_computer_time(computer_time, sync_times_computer, daq_sync) for computer_time in reward_times])

    coordinates = mat['fztrack']['tracks'][0][0][0].tolist()[0].T
    orientation = np.nan_to_num(mat['fztrack']['orientation'][0][0][0][0][0].flatten())
    angular_velocity = np.diff(orientation, prepend=0)

    prev_syllable = None
    syllable_bout_idx = 0

    nframes = len(syllable_by_frame)

    for frame in range(nframes):
        syll_id = syllable_by_frame[frame]

        curr_syllable = 'None' if syll_id not in syllable_index_to_name.keys() else syllable_index_to_name[syll_id]
        if curr_syllable != prev_syllable:
            if prev_syllable is not None:
                # finish processing the previous bout
                end_frame = frame-1
                bout_data["end_frame"] = end_frame
                bout_data["end_time"] = convert_cam_time(end_frame, sync_frames_foraging, daq_sync)
                bout_data["duration"] = bout_data["end_time"] - bout_data["start_time"]
                bout_data["end_x"] = float(coordinates[end_frame][0])
                bout_data["end_y"] = float(coordinates[end_frame][1])
                bout_data["next_syllable"] = curr_syllable
                trajectory_x = coordinates[bout_data["start_frame"]:bout_data["end_frame"]+1,0]
                trajectory_y = coordinates[bout_data["start_frame"]:bout_data["end_frame"]+1,0]
                bout_data["mean_orientation"] = np.mean(orientation[bout_data["start_frame"]:bout_data["end_frame"]+1])
                bout_data["mean_angular_velocity"] = np.sum(angular_velocity[bout_data["start_frame"]:bout_data["start_frame"]+int(bout_data["end_frame"]-bout_data["start_frame"]/2)])
                bout_data["mean_velocity"] = np.mean(np.sqrt(np.diff(trajectory_x)**2 + np.diff(trajectory_y)**2))
                bouts.append(bout_data)
                syllable_bout_idx += 1
                prev_syllable = bout_data["syllable"]

            bout_data = {}
            bout_data["syllable_bout_idx"] = syllable_bout_idx
            bout_data["start_frame"] = frame
            bout_data["start_time"] = convert_cam_time(frame, sync_frames_foraging, daq_sync)
            bout_data["syllable"] = curr_syllable
            bout_data["prev_syllable"] = prev_syllable
            prev_reward_ix = np.searchsorted(reward_times, bout_data["start_time"], side='right')-1
            try:
                bout_data["prev_reward"] = outcomes[prev_reward_ix]
            except:
                bout_data["prev_reward"] = None
            time_from_prev = bout_data["start_time"] - reward_times[prev_reward_ix]
            bout_data["time_from_prev_reward"] = None if time_from_prev < 0 else time_from_prev
            bout_data["time_to_reward"] = reward_times[prev_reward_ix+1] - bout_data["start_time"] if (prev_reward_ix < len(reward_times)-1) else None
            bout_data["next_reward"] = outcomes[prev_reward_ix+1] if (prev_reward_ix < len(reward_times)-1) else None
            bout_data["cumulative_food"] = total_food[prev_reward_ix]
            bout_data["cumulative_water"] = total_water[prev_reward_ix]
            bout_data["reward_index"] = prev_reward_ix
            bout_data["start_x"] = float(coordinates[frame][0])
            bout_data["start_y"] = float(coordinates[frame][1])

            prev_syllable = curr_syllable

    return pd.DataFrame(bouts)
