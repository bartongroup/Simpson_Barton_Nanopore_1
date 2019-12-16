import os

os.environ['KERAS_BACKEND'] = 'tensorflow'
from keras.models import load_model


MODELS = {
    'SQK-RNA001': 'sqkrna001_adapter_model.h5',
    'SQK-RNA002': 'sqkrna002_adapter_model.h5',
}


def load_adapter_model(model_fn=None, rna_kit='SQK-RNA001'):
    if model_fn is None:
        path = os.path.split(__file__)[0]
        model_fn = os.path.join(path,
                                'data',
                                MODELS[rna_kit])
    if not os.path.exists(model_fn):
        raise OSError('Model file {} does not exist'.format(model_fn))
    model = load_model(model_fn)
    return model
