import os
import gzip
import json
import logging

import numpy as np
import pandas as pd

import h5py
os.environ['KERAS_BACKEND'] = "tensorflow"
from keras import layers, models, optimizers, losses, callbacks
from keras.utils import Sequence
from keras import backend as K
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.utils import shuffle
from joblib import Parallel, delayed
import pysam
import click
import click_log

from .signal import get_noise_signal, random_dup_del, get_fast5_fiveprime
from .utils import get_fast5_read_id_mapping


EPOCHS = 100
BATCH_SIZE = 128
STEPS_PER_EPOCH = 100
VAL_STEPS_PER_EPOCH = 20

logger = logging.getLogger(__name__)
click_log.basic_config(logger)


class SquiggleSequence(Sequence):

    def __init__(self, pos_data, neg_data,
                 batch_size, steps_per_epoch,
                 random_frac=0.0625, squiggle_len=3000):
        self.pos_data = pos_data
        self.neg_data = neg_data
        self.batch_size = batch_size
        self.steps_per_epoch = steps_per_epoch
        self.random_frac = random_frac
        self.squiggle_len = squiggle_len
        self._augment_data()

    def _augment_data(self):
        n_pos = (self.batch_size * self.steps_per_epoch) // 2
        n_rand = int(n_pos * self.random_frac)
        n_neg = n_pos - n_rand
        pos_data_idx = np.random.randint(0, len(self.pos_data) - 1, size=n_pos)
        neg_data_idx = np.random.randint(0, len(self.pos_data) - 1, size=n_neg)
        X_data = np.concatenate(
            [self.pos_data[pos_data_idx],
             self.neg_data[neg_data_idx],
             [get_noise_signal(self.squiggle_len) for _ in range(n_rand)]]
        ).reshape(-1, self.squiggle_len, 1)
        y_data = np.concatenate(
            [[1] * n_pos, [0] * n_pos]
        ).reshape(-1, 1)
        #X_data_aug = np.array([
        #    random_dup_del(x) for x in X_data
        #])
        self.X_data_aug, self.y_data_aug = shuffle(X_data, y_data)

    def on_epoch_end(self):
        del self.X_data_aug
        del self.y_data_aug
        self._augment_data()

    def __len__(self):
        return self.steps_per_epoch

    def __getitem__(self, idx):
        idx = idx * self.batch_size
        X_batch = self.X_data_aug[idx: idx + self.batch_size]
        y_batch = self.y_data_aug[idx: idx + self.batch_size]
        return X_batch, y_batch
    

def generate_training_data(h5_file,
                           batch_size=128,
                           steps_per_epoch=25,
                           random_frac=0.0625,
                           test_split=0.1,
                           validation_split=0.1,
                           val_steps_per_epoch=5):
    '''
    Creates a training generator, a validation generator, and a test set from
    a HDF5 file containing pos_signals, neg_signals and neg_internal_signals.
    Data is augmented by addition of random noise and random duplicat
    '''
    logger.info('INFO: Reading in HDF data')
    with h5py.File(h5_file) as h5_file:
        pos = h5_file['pos_signals'][:]
        neg = np.concatenate([
            h5_file['neg_signals'][:],
            h5_file['neg_internal_signals'][:]
        ])
    logger.info('INFO: Splitting training and test datasets')
    squiggle_len = len(pos[0])
    pos_train, pos_test = train_test_split(pos, test_size=test_split)
    pos_train, pos_val = train_test_split(pos_train, test_size=validation_split)
    neg_train, neg_test = train_test_split(neg, test_size=test_split)
    neg_train, neg_val = train_test_split(neg_train, test_size=validation_split)

    # make test data 1:1
    if len(pos_test) < len(neg_test):
        neg_test = neg_test[:len(pos_test)]
    elif len(pos_test) > len(neg_test):
        pos_test = pos_test[:len(neg_test)]
    X_test = np.concatenate([pos_test, neg_test]).reshape(-1, squiggle_len, 1)
    y_test = np.concatenate([
        np.repeat(1, len(pos_test)),
        np.repeat(0, len(neg_test))
    ]).reshape(-1, 1)
    X_test, y_test = shuffle(X_test, y_test)
    logger.info('INFO: Building Sequence generators')
    train_sequence = SquiggleSequence(pos_train, neg_train,
                                      batch_size, steps_per_epoch,
                                      random_frac, squiggle_len)
    val_sequence = SquiggleSequence(pos_val, neg_val,
                                    batch_size, val_steps_per_epoch,
                                    random_frac, squiggle_len)

    return train_sequence, val_sequence, (X_test, y_test), squiggle_len


def conv_batch(prev, num_channels, kernel_size, name, conv_type=None):
    if conv_type is None:
        conv_type = layers.Conv1D
    prev = conv_type(num_channels,
                     kernel_size=kernel_size,
                     name=name + '_conv',
                     padding='same')(prev)
    prev = layers.BatchNormalization(name=name + '_batch_norm')(prev)
    return prev


def conv_batch_relu(prev, num_channels, kernel_size, name, conv_type=None):
    if conv_type is None:
        conv_type = layers.Conv1D
    prev = conv_batch(prev, num_channels, kernel_size, name, conv_type)
    prev = layers.Activation('relu', name=name + '_relu')(prev)
    return prev


def residual_block(prev, num_channels, name, conv_type=None, dim=1):
    shortcut = prev
    prev = conv_batch_relu(prev,
                           num_channels,
                           kernel_size=(5, ) * dim,
                           name=name + '_1',
                           conv_type=conv_type)
    prev = conv_batch(prev,
                      num_channels,
                      kernel_size=(5, ) * dim,
                      name=name + '_2',
                      conv_type=conv_type)

    shortcut = conv_batch(shortcut,
                          num_channels,
                          kernel_size=(1, ) * dim,
                          name=name + '_shortcut',
                          conv_type=conv_type)

    prev = layers.add([shortcut, prev], name=name + '_add_shortcut')
    prev = layers.Activation('relu', name=name + '_final_relu')(prev)
    return prev


def stack_convs(prev, num_channels, module_name, min_shape):
    convs = []
    i = 1
    while True:
        res_block_name = '{}_res_block_{}'.format(module_name, i)
        prev = residual_block(prev, num_channels, res_block_name, dim=1)
        curr_shape = tuple(prev.shape.as_list()[1:-1])
        if curr_shape <= min_shape:
            break
        convs.append(prev)
        prev = layers.MaxPool1D(2, name=res_block_name + '_max_pool', padding='same')(prev)
        i += 1
    return prev


def build_model(seq_length):
    input_layer = layers.Input(shape=(seq_length, 1))
    prev = layers.GaussianNoise(0.1)(input_layer)
    center = stack_convs(prev, 8, 'convs', min_shape=(24,))

    # binary output:
    prev = layers.Flatten()(center)
    prev = layers.Dense(16, activation='relu')(prev)
    prev = layers.Dropout(0.5)(prev)
    prev = layers.Dense(1, activation='sigmoid', name='binary_output')(prev)

    model = models.Model(input_layer, prev)
    model.compile(
        optimizer=optimizers.RMSprop(),
        loss=losses.binary_crossentropy,
    )
    return model


def train_model(model, train_gen, val_gen,
                epochs=25, verbose=2,
                steps_per_epoch=100,
                val_steps_per_epoch=20):
    lrd = callbacks.ReduceLROnPlateau(patience=3, min_delta=0.001)
    es = callbacks.EarlyStopping(patience=5)

    history = model.fit_generator(
        train_gen,
        steps_per_epoch=steps_per_epoch,
        epochs=epochs,
        verbose=verbose,
        validation_data=val_gen,
        validation_steps=val_steps_per_epoch,
        callbacks=[lrd, es],
        use_multiprocessing=False,
    )
    return model, history.history


def get_fiveprime_for_reads_in_bam(bam_fn, read_id_filemap, size=3000, internal=False, proc=28):
    with pysam.AlignmentFile(bam_fn) as bam:
        signals = Parallel(n_jobs=proc)(
            delayed(get_fast5_fiveprime)(
                r.query_name, read_id_filemap[r.query_name], size, internal)
            for r in bam.fetch()
        )
    if internal:
        read_ids, sig_lens, fiveprime_signals, internal_signals = zip(*signals)
        fiveprime_signals = [x for x in fiveprime_signals if x is not None and len(x) == size]
        internal_signals = [x for x in internal_signals if x is not None and len(x) == size]
        return np.asarray(fiveprime_signals), np.asarray(internal_signals)
    else:
        read_ids, sig_lens, fiveprime_signals = zip(*signals)
        fiveprime_signals = [x for x in fiveprime_signals if x is not None and len(x) == size]
        return np.asarray(fiveprime_signals)


def write_to_hdf5(h5_fn, pos_signals, neg_signals, neg_internal_signals):
    with h5py.File(h5_fn, 'w') as f:
        f.create_dataset('pos_signals', data=pos_signals)
        f.create_dataset('neg_signals', data=neg_signals)
        f.create_dataset('neg_internal_signals', data=neg_internal_signals)


def full_pass_training(h5_fn, squiggle_size, output_prefix, training_round):
    logger.info(f'INFO: Training pass {training_round}: Loading training data generators')
    train_gen, val_gen, test_data, _ = generate_training_data(
        h5_fn,
        batch_size=BATCH_SIZE,
        steps_per_epoch=STEPS_PER_EPOCH,
        val_steps_per_epoch=VAL_STEPS_PER_EPOCH
    )
    logger.info(f'INFO: Training pass {training_round}: Building new model')
    model = build_model(squiggle_size)
    logger.debug('DEBUG: Model Summary:')
    model.summary(print_fn=logger.debug)
    model, history = train_model(
        model, train_gen, val_gen,
        epochs=EPOCHS,
        steps_per_epoch=STEPS_PER_EPOCH,
        val_steps_per_epoch=VAL_STEPS_PER_EPOCH
    )
    X_test, y_test = test_data
    logger.info(f'INFO: Training pass {training_round} complete')
    y_pred = model.predict(X_test)
    logger.info('INFO: ROC AUC Score: {:.2f}'.format(roc_auc_score(y_test, y_pred)))
    history.update({
        'y_test': y_test.ravel().tolist(),
        'y_pred': y_pred.ravel().tolist()
    })
    # remove numpy types
    history = {k: [float(n) for n in v] for k, v in history.items()}
    model.save(f'{output_prefix}_model_wieghts_round{training_round}.h5')
    with gzip.open(
            f'{output_prefix}_training_history_round{training_round}.json',
            'wt', encoding="ascii") as f:
        json.dump(history, f)
    return model


@click.group()
def cli():
    pass


@cli.command()
@click_log.simple_verbosity_option(logger)
@click.option('-p', '--pos-bam-fn', required=True)
@click.option('-n', '--neg-bam-fn', required=True)
@click.option('-f', '--fast5-basedir', required=True)
@click.option('-s', '--summaries-basedir', required=True)
@click.option('-o', '--output-h5-fn', required=True)
@click.option('-l', '--squiggle-size', default=3000)
@click.option('--processes', default=1)
def build_training_data(pos_bam_fn, neg_bam_fn, fast5_basedir,
                        summaries_basedir, output_h5_fn,
                        squiggle_size, processes):
    read_id_filemap = get_fast5_read_id_mapping(
        summaries_basedir, fast5_basedir
    )
    logger.info(
        ('INFO: Found {} read_ids with associated Fast5'
        ' files in sequencing summaries').format(len(read_id_filemap)))
    pos_signals = get_fiveprime_for_reads_in_bam(
        pos_bam_fn, read_id_filemap, size=squiggle_size,
        internal=False, proc=processes
    )
    logger.info('INFO: Collected {} positive training examples'.format(len(pos_signals)))
    neg_signals, neg_signals_internal = get_fiveprime_for_reads_in_bam(
        neg_bam_fn, read_id_filemap, size=squiggle_size,
        internal=True, proc=processes
    )
    logger.info('INFO: Collected {} negative training examples'.format(len(neg_signals)))
    logger.info(
        'INFO: Collected {} negative internal training examples'.format(len(neg_signals_internal))
    )
    write_to_hdf5(
        output_h5_fn,
        pos_signals, neg_signals, neg_signals_internal
    )
    logger.info(
        f'INFO: Training data written to {output_h5_fn}'
    )


@cli.command()
@click_log.simple_verbosity_option(logger)
@click.option('-i', '--input-h5-fn', required=True)
@click.option('-o', '--output-prefix', required=True)
@click.option('-2/-1', '--two-pass-filt/--one-pass-filt', default=True)
@click.option('-l', '--squiggle-size', default=3000)
@click.option('--processes', default=1)
def run_training(input_h5_fn, output_prefix,
                 two_pass_filt, squiggle_size, processes):
    logger.info(
        'INFO: Running in {} pass mode'.format('two' if two_pass_filt else 'one')
    )

    model = full_pass_training(
        input_h5_fn,
        squiggle_size,
        output_prefix,
        training_round=1
    )
    if two_pass_filt:
        logger.info('INFO: Training pass 2')
        with h5py.File(f'{output_prefix}_training_data_round1.h5') as f:
            pos_signals = f['pos_signals'][:]
            neg_signals = f['neg_signals'][:]
            neg_signals_internal = f['neg_internal_signals'][:]

        preds = model.predict(neg_signals.reshape(-1, 3000, 1)).squeeze()
        neg_signals_filt = neg_signals[preds < 0.5]
        logger.info(
            'INFO: Filtered to {} negative training examples'.format(len(neg_signals_filt))
        )
        write_to_hdf5(
            f'{output_prefix}_training_data_round2.h5',
            pos_signals, neg_signals_filt, neg_signals_internal
        )
        logger.info(
            f'INFO: Training data written to {output_prefix}_training_data_round2.h5'
        )
        # retrain without false negatives
        model = full_pass_training(
            f'{output_prefix}_training_data_round2.h5',
            squiggle_size,
            output_prefix,
            training_round=2
        )

