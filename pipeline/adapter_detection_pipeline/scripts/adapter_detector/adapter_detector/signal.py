import os
import random

import noise
import numpy as np
from ont_fast5_api.multi_fast5 import MultiFast5File


def mad_scaling(signal):
    '''
    Scale a signal using median absolute deviation method
    '''
    shift = np.median(signal)
    scale = np.median(np.abs(signal - shift))
    return (signal - shift) / scale


def get_internal(signal, signal_size):
    sig_len = len(signal)
    if sig_len < signal_size:
        internal = np.zeros(size)
        internal[size - len(signal):] = signal
    elif sig_len < signal_size * 2:
        i = np.random.randint(0, sig_len - signal_size)
        internal = signal[i: i + signal_size]
    elif sig_len > 10000 + signal_size * 2:
        i = np.random.randint(10000, sig_len - signal_size * 2)
        internal = signal[i: i + signal_size]
    else:
        i = np.random.randint(0, sig_len - signal_size * 2)
        internal = signal[i: i + signal_size]
    return internal


def get_fast5_fiveprime(read_id, fast5_fns, signal_size, include_internal=False):
    '''
    Open fast5 file and return final signal_size measurements
    (corresponding to 5' end of an RNA signal). Signals
    are MAD scaled before the end is cropped.
    '''
    for fast5_fn in fast5_fns:
        with MultiFast5File(fast5_fn) as f5:
            try:
                read = f5.get_read(read_id)
            except KeyError as e:
                continue
            end = read.handle['Raw'].attrs['duration']
            signal = read.get_raw_data(scale=True, start=0, end=end)
            break
    else:
        return read_id, np.nan, np.empty(signal_size)
    signal = mad_scaling(signal)
    sig_len = len(signal)
    if sig_len >= signal_size:
        fiveprime = signal[sig_len - signal_size:]
    else:
        fiveprime = np.zeros(size)
        fiveprime[size - len(signal):] = signal
    if include_internal:
        internal = get_internal(signal, signal_size)
        return read_id, sig_len, fiveprime, internal
    return read_id, sig_len, fiveprime


def add_noise(mean, stdv):
    '''
    Add normally distributed noise to a signal (useful for data augmentation during training)
    '''
    mean = mean + np.random.normal(loc=0, scale=stdv, size=mean.shape)
    return mean


def get_noise_signal(signal_size, random_type=None):
    '''
    copied from DeepBinner balance module, creates one of four types of totally random
    negative signal for training.
    '''
    if random_type is None:
        random_type = random.choice(['flat', 'gaussian', 'multi_gaussian', 'perlin'])
    signal = []

    if random_type == 'flat':
        signal = [random.randint(0, 1000)] * signal_size

    elif random_type == 'gaussian':
        mean = random.uniform(300, 600)
        stdev = random.uniform(10, 500)
        signal = [int(random.gauss(mean, stdev)) for _ in range(signal_size)]

    elif random_type == 'multi_gaussian':
        while len(signal) < signal_size:
            length = random.randint(100, 500)
            mean = random.uniform(300, 600)
            stdev = random.uniform(10, 500)
            signal += [int(random.gauss(mean, stdev)) for _ in range(length)]

    elif random_type == 'perlin':
        octaves = random.randint(1, 4)
        step = random.uniform(0.001, 0.04)
        start = random.uniform(0.0, 1.0)
        factor = random.uniform(10, 300)
        mean = random.uniform(300, 600)
        signal = [int(mean + factor * noise.pnoise1(start + (i * step), octaves))
                  for i in range(signal_size)]

    signal = signal[:signal_size]
    return mad_scaling(add_noise(np.array(signal), 0.1))


def random_bincount_arr(size, arr_sum, replace=False):
    '''
    Creates an integer array with duplications and deletions
    which can be used to augment training data.
    '''
    while True:
        # something is dodgy with threading and np.random.choice
        idx = np.arange(int(size), dtype='i')
        idx = np.random.choice(idx, size=int(arr_sum), replace=replace).astype('i')
        if (idx >= 0).all():
            break
    return np.bincount(idx, minlength=int(size))


def random_dup_del(signal, min_frac=0, max_frac=0.33):
    '''
    Augment training data by duplicating some signal points and deleting others.
    Should help prevent overfitting of training data
    '''
    signal = signal.squeeze()
    sig_len = len(signal)
    frac = np.random.uniform(min_frac, max_frac)
    n_dups = sig_len // frac
    idx1 = random_bincount_arr(sig_len, n_dups, True)
    duplicated_signal = np.repeat(signal, idx1 + 1)
    idx2 = random_bincount_arr(sig_len + n_dups, sig_len, False).astype(bool)
    deleted_signal = duplicated_signal[idx2].copy()
    signal = deleted_signal.reshape(-1, 1)
    return signal

