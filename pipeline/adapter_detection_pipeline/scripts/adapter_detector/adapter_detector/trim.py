import os
import random
import numpy as np

os.environ['KERAS_BACKEND'] = 'tensorflow'
from keras.models import load_model

SEQ_OHE = {'A': [1, 0, 0, 0],
           'a': [1, 0, 0, 0],
           'C': [0, 1, 0, 0],
           'c': [0, 1, 0, 0],
           'G': [0, 0, 1, 0],
           'g': [0, 0, 1, 0],
           'T': [0, 0, 0, 1],
           't': [0, 0, 0, 1],
           'U': [0, 0, 0, 1],
           'u': [0, 0, 0, 1],
           'N': [0, 0, 0, 0]}


def zero_pad(seq):
    if len(seq) == 48:
        return seq
    rpad, r = divmod(48 - len(seq), 2)
    lpad = rpad + r
    return 'N' * lpad + seq + 'N' * rpad


def one_hot_sequence(seq, pad='zero'):
    ohe = []
    seq_len = len(seq)
    if seq_len < 48:
        seq = zero_pad(seq)
    elif seq_len > 48:
        raise ValueError('sequence is {} bases long'.format(seq_len))
    for base in seq:
        try:
            ohe.append(SEQ_OHE[base])
        except KeyError:
            ohe.append(SEQ_OHE['N'])
    return np.array(ohe)


class Trimmer():
    '''
    Use a model trained on sequence level to identify the length of the
    input sequence (the first 48nt of the basecalled RNA) which is adapter,
    and trim it off
    '''

    def __init__(self, model_fn=None):
        if model_fn is None:
            path = os.path.split(__file__)[0]
            model_fn = os.path.join(path,
                                    'data',
                                    'trimmer_model.h5')
        self.model = load_model(model_fn)


    def trim_adapters(self, read_batch, needs_trimming):
        seq_ohe = np.array([one_hot_sequence(read.sequence[:48]) for read in read_batch])
        idx = self.model.predict(seq_ohe).ravel()
        idx = np.round(idx).astype(np.int)
        idx = np.maximum(idx, 0)
        # leave ambiguous sequences in place
        idx[idx > 48] = 0
        for i, read, t in zip(idx, read_batch, needs_trimming):
            if t:
                read.sequence = read.sequence[i:]
                read.quality = read.quality[i:]
        return read_batch