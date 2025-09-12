import pandas as pd
import numpy as np
from scipy.stats import binned_statistic
import os
import pickle
import glob
from joblib import Parallel, delayed
from matplotlib.colors import ListedColormap

from anatomy import *
from lfp import *

## Set this to histology export path for convenience ##
## Alternatively pass a path to "load_histology()" ##
HISTOLOGY_PATH = "Z:/motiv/data/histology/hist_export_buridan_alt"

class MultiprobeEphysExperiment(object):
    """
    Object representing a combined recording with behavior, histology, and spiking data for a single recording session.
    
    This object can contain multiple simultaneously recorded sorted datasets + multiple 
    behavioral sessions. See jupyter notebook for example usage.

    For convenience, this class may be sub-classed to add experiment-specific
    behavioral data or event processing, such as trials and other event triggers.

    """
    def __init__(self, penetration=None, recording=None, mouse_name=None):
        self._recording = None
        self._behavior = [] # list of SessionBehavior objects containing behavior information for each session
        self._video = [] # VideoData objects with PCs representing animal motion during behavior, for each session.
        self._histology = None
        self._unit_data = pd.DataFrame()
        self._units = None
        self._nrecording = 0

        # database values
        self._penetration = penetration
        self._recording = recording
        self._mouse_name = mouse_name
        if penetration is not None:
            self._tract_names = penetration.split(';')
        else:
            self._tract_names = []

        self._unit_locs = None
        self._unit_areas = None
        self._unit_area_ids = None

        self._atlas = None

    def load_histology(self, atlas, histology_path=None):

        if histology_path is None:
            histology_path = HISTOLOGY_PATH
        
        fname = os.path.join(histology_path,
                    self._mouse_name + "_" + self._recording + ".npz")
        if os.path.exists(fname):
            S = np.load(fname)
            self._atlas = atlas
            self._unit_locs = np.squeeze(S["transformed_locs"])
            self._chan_locs = np.squeeze(S["transformed_chans"])
            self._chan_probe_id = np.squeeze(S["chan_probe_id"])
        else:
            print('Path not found: %s' % fname )

        # Handle case where not every ephys probe has anatomy
        # by setting locations to NaN for that probe.
        # These units should be filtered out for any anatomical analysis.
        if self._probe_id is not None:
            ephys_probes = np.unique(self._probe_id)
            anatomy_probes = np.unique(self._chan_probe_id)
            missing_probe_anatomy = np.sort([p for p in ephys_probes if p not in anatomy_probes])
            if len(missing_probe_anatomy) != 0:
                for probe in missing_probe_anatomy:
                    chan_probe_id = self._probe_id[self._probe_id==probe]
                    unit_locs = np.full((len(chan_probe_id),3), np.nan)
                    chan_locs = np.full((len(chan_probe_id),3), np.nan)
                    self._chan_probe_id = np.hstack([self._chan_probe_id, chan_probe_id])
                    self._unit_locs = np.vstack([self._unit_locs, unit_locs])
                    self._chan_locs = np.vstack([self._chan_locs, chan_locs])

    def load_units_from_ks(self, fpath, verbose=True, load_lfp=False):
        """
        Load units from numpy array, then compute firing rates.
        NOTE: The KS data should have been aligned to sync events, so we can
        assume that everything is in the same timeframe.
        """
        if ".npz" in fpath:
            data = np.load(fpath, allow_pickle=True)
        else:
            with open(fpath, 'rb') as f:
                data = pickle.load(f)
        
        self._spike_times = data._spike_times
        self._clusts = data._clusts
        self._clust_id = data._clust_id
        self._clust_depths = data._clust_depths
        self._probe_id = data._probe_id
        self._metrics = data._metrics

        if load_lfp:
            self._lfp = data._lfp
            self._lfp_fs = data._lfp_fs
            self._lfp_chans = data._lfp_chans
            self._lfp_probe_ids = data._lfp_probe_ids
        else:
            self._lfp = None
            self._lfp_fs = None
            self._lfp_chans = None
            self._lfp_probe_ids = None

        # convenience data for performant spike extraction
        self.__clusts_argsorted = np.argsort(self._clusts)
        self.__clusts_sorted = self._clusts[self.__clusts_argsorted]
        # event sync data
        self._events = data._events

        self._ntot = data._ncell
        self._max_t = self._spike_times.max()
        self.get_good_units()

    def load_lfp_raw(self, raw_dir, processed_dir, spatial_downsample=4, temporal_downsample=10, do_save=True, use_cache=True):

        import scipy
        import sglx_util
        import lfp

        probes = np.unique(self._probe_id)

        all_lfp= []
        all_fs = []
        lfp_chans = []
        lfp_probe_ids = []

        self._sync_times = {}

        save_path = glob.glob(os.path.join(processed_dir, self._mouse_name, f"*{self._recording}*", "*", f"*imec{probes[0]}*") + f"*chanMap.mat")[0]
        save_path = os.path.abspath(os.path.join(os.path.dirname(save_path),"..","lfp_downsampled.pkl"))

        if use_cache and os.path.exists(save_path):
            self.load_lfp_pkl(save_path)
            return

        for probe in probes:
            print(probe)
            raw_path = os.path.join(raw_dir, self._mouse_name, f"*{self._recording}*", "*", f"*imec{probe}*")
            processed_path = os.path.join(processed_dir, self._mouse_name, f"*{self._recording}*", "*", f"*imec{probe}*")
            
            lf_fn = glob.glob(raw_path + f"*lf.bin")[0]
            chanmap_fn = glob.glob(processed_path + "*chanMap.mat")[0]
            connected = np.squeeze(scipy.io.loadmat(chanmap_fn)["connected"])
            lfp_data = sglx_util.Reader(lf_fn)
            fs = lfp_data.fs
            if fs == None:
                fs = 2500.
            # subsample LFP data spatially and temporally -- take every 4th channel, every 10 samples
            lfp_sub = lfp.subsample_lfp(lfp_data._raw, np.arange(384)[connected.astype('bool')][::spatial_downsample], temporal_downsample)
            chan_ix = scipy.io.loadmat(chanmap_fn)['chanMap0ind'].flatten()[connected.astype('bool')].astype('int')[::spatial_downsample]
            lfp_data.close()

            # Load per-probe sync times
            sync_fn = glob.glob(processed_path + "*ap.SY*")[0]
            # Load each line of sync_fn into a numpy array
            sync_times = np.loadtxt(sync_fn)
            self._sync_times[probe] = sync_times

            all_lfp.append(lfp_sub)
            lfp_chans.append(chan_ix)
            lfp_probe_ids.append(probe*np.ones(len(chan_ix)))

            all_fs.append(fs/temporal_downsample)

        self._lfp = all_lfp
        self._all_fs = all_fs
        self._lfp_chans = np.hstack(lfp_chans)
        self._lfp_probe_ids = np.hstack(lfp_probe_ids)

        import pickle

        if do_save:
            with open(save_path, 'wb') as f:
                pickle.dump({
                    'lfp': self._lfp, 
                    'lfp_fs': self._all_fs, 
                    'lfp_chans': self._lfp_chans, 
                    'lfp_probe_ids': self._lfp_probe_ids, 
                    'sync_times': self._sync_times
                }, f, protocol=pickle.HIGHEST_PROTOCOL)
        #np.savez(save_path, lfp=self._lfp, lfp_fs=fs, lfp_chans=self._lfp_chans, lfp_probe_ids=self._lfp_probe_ids, sync_times=self._sync_times, allow_pickle=True)

    def load_lfp_pkl(self, fpath):
        # load pickled LFP data
        with open(fpath, 'rb') as f:
            data = pickle.load(f)
        self._lfp = data["lfp"]
        self._lfp_fs = data["lfp_fs"]
        self._lfp_chans = data["lfp_chans"]
        self._lfp_probe_ids = data["lfp_probe_ids"]
        self._sync_times = data["sync_times"]

    def get_good_units(self, isi_viol=0.1, snr=1.5, nspikes=500, noise=0, qc='bombcell'):
        m = self._metrics
        if qc=='bombcell':
            good_annot = m["good"].astype('bool')
        else:
            good_annot = m["noise"]==noise
        #mask = (m["isi_viol"]<isi_viol) & (m["snr"]>=1.5) & (m["n_spikes"]>nspikes) & good_annot
        mask = (m["snr"]>=1.5) & (m["n_spikes"]>nspikes) & good_annot
        self._good = np.unique(m[mask]["cluster_id"])
        self._ncell = len(self._good)
        self._good_mask = mask

    def get_expt_name(self):
        return self._penetration + "_" + self._mouse_name + "_" + self._recording

    def get_ncells(self):
        return self._ncell
 
    def get_clust_ids(self):
        return self._clust_id

    def get_clust_groups(self):
        """
        Return cluster groups: good, MUA, noise
        """
        return self._clust_groups

    def get_maxt(self):
        return self._max_t

    def get_lfp_for_probe(self,probe_id):
        return self._lfp[:, self._lfp_probe_ids==probe_id]

    #
    # METHODS RELATED TO ANATOMY
    #

    def get_brain_area_map(self, use_higher=False, locs=None):
        mask = True if locs is None else False
        ids = self.get_brain_areas(info='id', mask=mask, locs=locs)
        acronyms = self.get_brain_areas(info='acronym', mask=mask, locs=locs)
        name_to_id = {}
        id_to_name = {}
        for i,a in enumerate(acronyms):
            name_to_id[a] = ids[i]
            id_to_name[ids[i]] = a
        return name_to_id, id_to_name

    def get_brain_areas(self, info='id', level='bottom', mask=True, locs=None):
        if locs is None:
            locs = self._unit_locs
        # NOTE THAT THESE ARE SORTED BY DEPTH IN REVERSE ORDER!!
        # Split locs in to those with nans and those without.
        # For locs with nan values, return NA for acronym and -1 for id
        nan_locs = locs[np.where(np.any(np.isnan(locs),1))]
        nonan_locs = locs[np.where(np.all(~np.isnan(locs),1))]
        areas = self._atlas.get_info_for_points(nonan_locs, element="id", increase_level=False)
        areas = np.array([self._atlas.get_acronym(a, level=level) for a in areas])
        areas = np.hstack([areas, np.full(len(nan_locs), 'NA')])

        acronym_map = self._atlas._tree.get_id_acronym_map()

        if info == 'id':
            # doing this a second time to get the area id at the appropriate "level"
            areas = np.array([-1 if ar=='NA' else acronym_map[ar] for ar in areas])
        if mask:
            return areas[self._good_mask]
        else:
            return areas

    def get_areas_sorted_hierarchically(self, locs=None):
        """ Get index of good units sorted hierarchically by brain region subdivision
        """
        mask = True if locs is None else False
        area_ids = self.get_brain_areas(info='id', level='bottom', mask=mask, locs=locs)
        area_ix = np.arange(len(area_ids))
        areas_top = self.get_brain_areas(info='acronym', level='top', mask=mask, locs=locs)
        area_hierarchy_sorted = np.hstack([area_ix[areas_top==region][np.argsort(area_ids[areas_top==region])] for region in np.unique(areas_top)])
        return area_hierarchy_sorted

    def get_area_colors(self, area_ids):
        """
        Get Allen CCF colors of areas by id, with unassigned areas mapped to white.

        Example usage: plt.plot(regional_psth, color=get_area_colors([region_id]))
        """
        area_colors = np.array([(255,255,255) if aid==-1 else self._atlas._tree.get_colormap()[aid] for aid in area_ids])/255
        return area_colors

    def get_area_colorbar(self, area_ids):
        """
        Given a list of area ids, generate a colorbar (for plotting by imshow)
        and an associated colormap. Uses colors from the Allen CCF.

        Example usage: plt.imshow(colorbar, cmap=cm, aspect='auto')
        Example usage: plt.imshow(colorbar[sorted_idx,:], cmap=cm, aspect='auto')
        """
        id_color = []
        distinct_color = []
        for i,aid in enumerate(area_ids):
            if aid not in id_color:
                id_color.append(aid)
                if aid == -1:
                    distinct_color.append(np.array([1.,1.,1.]))
                else:
                    distinct_color.append(np.array(self._atlas._tree.get_colormap()[aid])/255)
        id_color = dict(zip(id_color, range(len(id_color))))

        colors = np.pad(np.array(distinct_color), [0,1], constant_values=1)[:-1,:]
        cm = ListedColormap(colors)
        colorbar = np.expand_dims([id_color[aid] for aid in area_ids],1)
        return colorbar, cm

    def preprocess_lfp(self, do_notch=False):
        """ This function preprpocesses the per-probe LFP data in the following ways:
        1. Remove "NA" channels (which should include those that are out of the brain)
        2. Common median reference
        3. Zscore each channel across time
        4. Notch filter to remove 60hz and 120hz noise
        """
        from scipy.stats import zscore

        def notch(sig, fs, notch_freq, q_fac=30):
            from scipy import signal as ss
            '''
            sig (np.ndarray): time x channels

            TODO: eventually add harmonics when we use >200 hz sampled data
            '''
            b_notch, a_notch = ss.iirnotch(notch_freq, q_fac, fs)
            return ss.filtfilt(b_notch, a_notch, sig, axis=0)
        
        # Restrict probes to only the anatomically reconstructed ones
        # Since we don't want LFP data from unused probes sitting in air
        valid_probes = [p for p in np.unique(self._chan_probe_id) if not np.alltrue(np.isnan(self._chan_locs[self._chan_probe_id==p]))]

        lfp_areas = self.get_brain_areas(mask=False, info='acronym', locs=self.get_locs_for_lfp())

        lfp_p = {}
        lfp_p_areas = {}
        probes = np.unique(self._lfp_probe_ids)

        for i,probe in enumerate(probes):
            if probe in valid_probes:
                probe_mask = np.array([pid for pid in self._lfp_probe_ids if pid in valid_probes])==probe
                area_mask = lfp_areas != 'NA'
                channel_mask = area_mask[probe_mask]
                car = self._lfp[i][:,channel_mask] - np.expand_dims(np.median(self._lfp[i][:,channel_mask],1),1)
                car = zscore(car, 0)
                if do_notch:
                    car = notch(car, 250., 60, q_fac=30)
                    car = notch(car, 250., 120, q_fac=30)
                lfp_p[probe] = car
                lfp_p_areas[probe] = lfp_areas[np.logical_and(area_mask, probe_mask)]

        self._lfp_p = lfp_p
        self._lfp_p_areas = lfp_p_areas
        return

    def get_locs_for_lfp(self):
        valid_probes = [p for p in np.unique(self._chan_probe_id) if not np.alltrue(np.isnan(self._chan_locs[self._chan_probe_id==p]))]
        return np.concatenate([self._chan_locs[self._chan_probe_id==p][self._lfp_chans[self._lfp_probe_ids==p].astype('int')] for p in np.unique(self._chan_probe_id) if p in valid_probes])

    def adjust_probe_time(self, t, probe):
        """Adjust an event time (on the reference clock) to the time on the target probe's clock.
        t: event time on reference
        p: probe id
        """
        # sync pulse times on reference clock
        ref = self._sync_times[0]
        # sync pulse times on target clock (must match number of ref sync pulses)
        target = self._sync_times[probe]

        prev_sync_ix = np.searchsorted(ref,t)
        offset = t-ref[prev_sync_ix]
        return target[prev_sync_ix]+offset

    def lfp_per_event(self, event_times, fwd_t, rev_t):
        """ Return LFP data for each event, for each probe, windowed around the event
        event_times: time of each event
        fwd_t: forward time in seconds
        rev_t: reverse time in seconds
        """
        dat = []
        all_lfp_areas = []
        for pix, probe in enumerate(np.array(list(self._lfp_p.keys())).astype('int')):
            probe_times = np.arange(len(self._lfp_p[probe]))/self._lfp_fs[pix]
            event_dat = []
            for e in event_times:
                e = self.adjust_probe_time(e, probe)
                start = np.searchsorted(probe_times, e-rev_t)
                total_time = int((fwd_t+rev_t)*self._lfp_fs[pix])
                stop = start+total_time
                event_dat.append(self._lfp_p[probe][start:stop,:].T)
            event_dat = np.swapaxes(np.stack(event_dat),0,1)
            dat.append(event_dat)
            all_lfp_areas.append(self._lfp_p_areas[probe])
        dat = np.concatenate(dat, 0)
        all_lfp_areas = np.concatenate(all_lfp_areas)
        return dat, all_lfp_areas
    
    #    
    # METHODS RELATING TO GETTING SPIKES
    # 
    def firing_rates_over_recording(self, bin_size=0.01, smooth_size=None):
        """
        Return a matrix containing the binned firing rate for each neuron 
        over the whole recording session.
        Returns the matrix and an array of the start time of each bin of 
        the matrix.             
        :input bin_size: bin size in seconds
        """
        nbins = self._max_t / float(bin_size)
        t = np.arange(0, self._max_t, bin_size)
        fr = np.zeros((len(self._good), len(t) - 1)) 
        for idx, clu in enumerate(self._good):
            curr_spikes = self.get_spikes_for_clust(clu).flatten()
            # array of ones of the same size as the spike times, 
            # to compute the count
            ones = np.ones_like(curr_spikes) 
            curr_fr, bin_edges, _ = binned_statistic(curr_spikes, ones, 
                statistic="sum", bins=nbins, range=(0, self._max_t))
            fr[idx, :] = curr_fr / bin_size
        if smooth_size:
            fr = moving_avg_filter_2d(fr, smooth_size)

        return fr, bin_edges

    def get_spikes_per_event(self, event_times,fwd_t, rev_t=0, shift=True, use_both=False):
        """
        Return list of nCell of list of nEvent of spikes per event
        Args:
            event_times: time in seconds of each event to obtain spikes from around
            fwd_t (float): length of event in seconds
            rev_t (float, optional): preceding time in seconds
            shift (bool, optional): whether to shift time of each spike relative to the event time. Default is True. NOTE: Corrected for rev_t so min spike time is 0, max is rev_t+fwd_t
        Returns:
            spikes_per_cell: list of spikes per event_times, per cell (Ncell list of Ntrials)
        """
        if use_both:
            cell_ids = np.hstack((self._good, self._mua))
        else:
            cell_ids = self._good
        spikes_per_cell = [None for i in range(len(cell_ids))]
        for i in range(len(spikes_per_cell)):
            spikes_per_cell[i] = [None]*len(event_times)
        # only operate on "good" cells for now (ignore MUA)
        for i,idx in enumerate(cell_ids):
            spikes_for_clust = np.sort(self.get_spikes_for_clust(idx))
            for t,s in enumerate(event_times):
                start,stop = np.searchsorted(spikes_for_clust,[s-rev_t, s+fwd_t])
                spikes_per_cell[i][t] = spikes_for_clust[start:stop]
                #spikes_per_cell[i][t] = spikes_for_clust[np.logical_and(spikes_for_clust >= (s - rev_t), spikes_for_clust < (s + fwd_t))] 
                if shift:
                    spikes_per_cell[i][t] -= s
                    spikes_per_cell[i][t] += rev_t # make spikes start at 0
        return spikes_per_cell
 
    # METHODS RELATED TO SPATIAL LOCATION ON PROBE
    #
    def get_clusters_sorted_by_depth(self, clust_type="good"):
        """
        Computes good units sorted by depth along electrode.
        Returns:
             cluster_id: sorted list of ids of clusters
             sorted_idx: indices of each cluster (e.g. to rearrange firing rate matrix)
        """
        if clust_type == "both":
            sorted_idx = np.argsort(self._clust_depths)
        elif clust_type == "mua":
            raise NotImplementedError
            #sorted_idx = np.argsort(self._clust_depths[self._clust_groups==1])
        elif clust_type == "good":
            sorted_idx = np.argsort(self._clust_depths[self._good])
        return self._clust_id[sorted_idx], sorted_idx
    
    def depth_sort_probe(self, probe):
        return self._clust_depths[(data._good) & (data._probe_id==probe)].argsort()

    def get_cluster_depths(self, clust_type="good"):
        if clust_type == "both":
            return self._clust_depths
        elif clust_type == "mua":
            raise NotImplementedError
            #return self._clust_depths[self._clust_groups==1]
        elif clust_type == "good":
            return self._clust_depths[self._good]

    #
    def get_spikes_for_clust(self,clust_id,use_timestamps=False):
        """
        Return timestamps for spikes from a particular cluster, in s
        """
        b, e = np.searchsorted(self.__clusts_sorted, [clust_id, clust_id+1])
        return self._spike_times[self.__clusts_argsorted[b:e]]
        #return self._spike_times[self._clusts == clust_id]


    def get_nearest_event(self, events, t):
        return events[np.argmin(np.abs(events-t))]

# Utilities

def moving_avg_filter(t,win_size=10,causal=True):
    if causal:
        filtered = np.append(t[:win_size-1], pd.Series(t).rolling(window=win_size).mean().iloc[win_size-1:].values)
    else:
        filtered = np.append(pd.Series(t).rolling(window=win_size).mean().iloc[win_size-1:].values, t[-win_size+1:])
    assert t.size==filtered.size
    return filtered

def moving_avg_filter_2d(fr, win_size=10, causal=True):
    assert len(fr.shape)==2
    return np.vstack(Parallel(n_jobs=-1)(delayed(moving_avg_filter)(fr[i,:], win_size=win_size, causal=causal) for i in np.arange(fr.shape[0])))

def moving_avg_filter_3d(fr, win_size=10, causal=True):
    assert len(fr.shape)==3
    return np.stack(Parallel(n_jobs=-1)(delayed(moving_avg_filter_2d)(fr[:,i,:], win_size=win_size, causal=causal) for i in np.arange(fr.shape[1])),1)

def find(x):
    return np.argwhere(x).flatten()