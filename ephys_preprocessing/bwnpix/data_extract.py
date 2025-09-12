import sys
import os
from load import *
from multiprobe import *
from bwnpix.data_loading import *
import scipy
from pylab import *
import pandas as pd
import scipy.io
sys.path.append(os.path.join(os.path.dirname(__file__), '../brainwide-npix/bwnpix'))

# Define paths for data storage and export
# Define paths for data storage and retrieval
cache_root = '/path/to/cache_root/'
behavior_path = '/path/to/npix_behavior/'
fig_path = '/path/to/figures/'
data_export_path = '/path/to/export/'
ephys_path = '/path/to/ephys/localProcessedData/'


def extract_recording_data(all_data):
    """
    Extracts recording data and saves each variable as a MATLAB .mat file.

    Args:
        all_data (list): A list of data objects containing recording information.
    """
    for data in all_data:
        recording_name = data._recording
        recording_fr, recording_bin_size = data.firing_rates_over_recording(bin_size=0.01, smooth_size=None)
        recording_info = {
            'recording_info': {
                'mouse_name': data._mouse_name,
                'recording_name': data._recording,
                'beh_comp_sync': data._events,
                'spike_times': data._spike_times,
                'firing_rates': recording_fr,
                'fr_bin_starts': recording_bin_size
            }
        }
        # Save each variable in recording_info as a separate .mat file
        for key, value in recording_info.items():
            file_path = os.path.join(data_export_path, f'{recording_name}_recording_info.mat')
            scipy.io.savemat(file_path, {key: value})

def spike_times_by_neuron_npy(all_data):
    """
    Save the spike times for all good units from a recording as .npy files.
    
    Args:
        all_data (list): A list of data objects containing recording information.
    """
    for data in all_data:
        recording_name = data._recording
        numgood = data._good.size
        spikes_neuron = np.empty((numgood,), dtype=object)
        for idx, clu in enumerate(data._good):
            curr_spikes = np.sort(data.get_spikes_for_clust(clu).flatten())
            spikes_neuron[idx] = curr_spikes
        
        # Save spikes_neuron as a .npy file
        file_path = os.path.join(data_export_path, f'{recording_name}_spikes_neuron.npy')
        print(f'{recording_name} is saving')
        np.save(file_path, spikes_neuron)

def spike_times_by_neuron_mat(all_data):
    """
    Save the spike times for all good units from a recording as .mat files.
    
    Args:
        all_data (list): A list of data objects containing recording information.
    """
    for data in all_data:
        recording_name = data._recording
        numgood = data._good.size
        spikes_neuron = np.empty((numgood,), dtype=object)
        for idx, clu in enumerate(data._good):
            curr_spikes = np.sort(data.get_spikes_for_clust(clu).flatten())
            spikes_neuron[idx] = curr_spikes
        
        # Save spikes_neuron as a .mat file
        file_path = os.path.join(data_export_path, f'{recording_name}_spikes_neuron.mat')
        print(f'{recording_name} is saving')
        scipy.io.savemat(file_path, {'spikes_neuron': spikes_neuron})

def extract_regions_per_session(all_data, recording_names):
    """
    Extracts region information per session and saves it as CSV files.

    Args:
        all_data (list): A list of data objects containing session information.
        recording_names (list): A list of recording names corresponding to the data objects.
    """
    for d, data in enumerate(all_data):  # Loop through all data entries
        locs = data._unit_locs[data._good_mask]  # Get unit locations for good units
        recorded_areas = data.get_brain_areas(info='id')  # Get recorded brain area IDs
        region_names = data.get_brain_areas(info='acronym')  # Get region acronyms
        all_region_cbar, all_recorded_cm = data.get_area_colorbar(recorded_areas)
        
        # Convert lists to DataFrames
        df_region_names = pd.DataFrame(region_names)
        df_recorded_areas = pd.DataFrame(recorded_areas)
        df_locs = pd.DataFrame(locs.tolist())
        df_all_region_cbar = pd.DataFrame(all_region_cbar)
        df_all_recorded_cm = pd.DataFrame(all_recorded_cm.colors)

        # Export DataFrames to CSV files for each recording name
        recording_name = recording_names[d]
        df_region_names.to_csv(os.path.join(data_export_path, f'{recording_name}_region_names.csv'), index=False)
        df_recorded_areas.to_csv(os.path.join(data_export_path, f'{recording_name}_recorded_areas.csv'), index=False)
        df_locs.to_csv(os.path.join(data_export_path, f'{recording_name}_locs.csv'), index=False)
        df_all_region_cbar.to_csv(os.path.join(data_export_path, f'{recording_name}_all_region_cbar.csv'), index=False)
        df_all_recorded_cm.to_csv(os.path.join(data_export_path, f'{recording_name}_all_recorded_cm.csv'), index=False)
    return

def extract_metrics_per_session(all_data, recording_names):
    """
    Extracts metrics information per session and saves it as CSV files.

    Args:
        all_data (list): A list of data objects containing session information.
        recording_names (list): A list of recording names corresponding to the data objects.
    """
    for d, data in enumerate(all_data):  # Loop through all data entries
        metrics = data._metrics  # Get metrics information
        metrics['good_2_use'] = data._good_mask.astype(int)  # Add good_2_use logical vector (1 for True, 0 for False)
        # Convert metrics to DataFrame
        df_metrics = pd.DataFrame(metrics)
        # Export DataFrame to CSV file 
        recording_name = recording_names[d]
        df_metrics.to_csv(os.path.join(data_export_path, f'{recording_name}_metrics.csv'), index=False)
    return

def printlocs(all_data):
    for data in all_data:
        recorded_areas = data.get_brain_areas(info='id')  # Get recorded brain area IDs
        indices = [index for index, value in enumerate(recorded_areas) if value == -1]
        print(f'Indices with id -1: {indices}')
        for idx in indices:
            print(f'Unit locs for index {idx}: {data._unit_locs[idx]}')
