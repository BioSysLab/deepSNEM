from __future__ import division, print_function
import numpy as np
from numpy import inf, ndarray
import pandas as pd
import tensorflow as tf
import os
import random
import keras
import sklearn
from sklearn.model_selection import KFold
from sklearn import metrics
import re
from keras import optimizers
from keras import losses
from keras import regularizers
import keras.backend as K
from keras.models import model_from_json
from keras.models import load_model
from tempfile import TemporaryFile
from keras import layers
from keras.layers import Input, Dense, merge, BatchNormalization, Activation
from keras.layers import deserialize as layer_from_config
from rdkit import Chem
from functools import partial
from multiprocessing import cpu_count, Pool
from keras.utils.generic_utils import Progbar
from copy import deepcopy
from NGF.utils import filter_func_args, mol_shapes_to_dims
import NGF.utils


def temporal_padding(x, padding=(1, 1), padvalue = 0):
  """Pads the middle dimension of a 3D tensor.
  Arguments:
      x: Tensor or variable.
      padding: Tuple of 2 integers, how many zeros to
          add at the start and end of dim 1.
  Returns:
      A padded 3D tensor.
  """
  assert len(padding) == 2
  pattern = [[0, 0], [padding[0], padding[1]], [0, 0]]
  return tf.pad(x, pattern,constant_values = padvalue)
  
  
  
  
def neighbour_lookup(atoms, edges, maskvalue=0, include_self=False):
    ''' Looks up the features of an all atoms neighbours, for a batch of molecules.

    # Arguments:
        atoms (K.tensor): of shape (batch_n, max_atoms, num_atom_features)
        edges (K.tensor): of shape (batch_n, max_atoms, max_degree) with neighbour
            indices and -1 as padding value
        maskvalue (numerical): the maskingvalue that should be used for empty atoms
            or atoms that have no neighbours (does not affect the input maskvalue
            which should always be -1!)
        include_self (bool): if True, the featurevector of each atom will be added
            to the list feature vectors of its neighbours

    # Returns:
        neigbour_features (K.tensor): of shape (batch_n, max_atoms(+1), max_degree,
            num_atom_features) depending on the value of include_self

    # Todo:
        - make this function compatible with Tensorflow, it should be quite trivial
            because there is an equivalent of `T.arange` in tensorflow.
    '''

    # The lookup masking trick: We add 1 to all indices, converting the
    #   masking value of -1 to a valid 0 index.
    masked_edges = edges + 1
    # We then add a padding vector at index 0 by padding to the left of the
    #   lookup matrix with the value that the new mask should get
    masked_atoms = temporal_padding(atoms, (1,0), padvalue=maskvalue)


    # Import dimensions
    atoms_shape = K.shape(masked_atoms)
    batch_n = atoms_shape[0]
    lookup_size = atoms_shape[1]
    num_atom_features = atoms_shape[2]

    edges_shape = K.shape(masked_edges)
    max_atoms = edges_shape[1]
    max_degree = edges_shape[2]

    # create broadcastable offset
    offset_shape = (batch_n, 1, 1)
    offset = K.reshape(tf.keras.backend.arange(stop=batch_n,start=0, dtype='int32'), offset_shape)
    offset *= lookup_size

    # apply offset to account for the fact that after reshape, all individual
    #   batch_n indices will be combined into a single big index
    flattened_atoms = K.reshape(masked_atoms, (-1, num_atom_features))
    flattened_edges = K.reshape(masked_edges + offset, (batch_n, -1))

    # Gather flattened
    flattened_result = K.gather(flattened_atoms, flattened_edges)

    # Unflatten result
    output_shape = (batch_n, max_atoms, max_degree, num_atom_features)
    output = K.reshape(flattened_result, output_shape)

    if include_self:
        return K.concatenate([tf.expand_dims(atoms, axis=2), output], axis=2)
    return output
	
	
	
class NeuralGraphHidden(layers.Layer):
    ''' Hidden Convolutional layer in a Neural Graph (as in Duvenaud et. al.,
    2015). This layer takes a graph as an input. The graph is represented as by
    three tensors.

    - The atoms tensor represents the features of the nodes.
    - The bonds tensor represents the features of the edges.
    - The edges tensor represents the connectivity (which atoms are connected to
        which)

    It returns the convolved features tensor, which is very similar to the atoms
    tensor. Instead of each node being represented by a num_atom_features-sized
    vector, each node now is represented by a convolved feature vector of size
    conv_width.

    # Example
        Define the input:
        ```python
            atoms0 = Input(name='atom_inputs', shape=(max_atoms, num_atom_features))
            bonds = Input(name='bond_inputs', shape=(max_atoms, max_degree, num_bond_features))
            edges = Input(name='edge_inputs', shape=(max_atoms, max_degree), dtype='int32')
        ```

        The `NeuralGraphHidden` can be initialised in three ways:
        1. Using an integer `conv_width` and possible kwags (`Dense` layer is used)
            ```python
            atoms1 = NeuralGraphHidden(conv_width, activation='relu', bias=False)([atoms0, bonds, edges])
            ```
        2. Using an initialised `Dense` layer
            ```python
            atoms1 = NeuralGraphHidden(Dense(conv_width, activation='relu', bias=False))([atoms0, bonds, edges])
            ```
        3. Using a function that returns an initialised `Dense` layer
            ```python
            atoms1 = NeuralGraphHidden(lambda: Dense(conv_width, activation='relu', bias=False))([atoms0, bonds, edges])
            ```

        Use `NeuralGraphOutput` to convert atom layer to fingerprint

    # Arguments
        inner_layer_arg: Either:
            1. an int defining the `conv_width`, with optional kwargs for the
                inner Dense layer
            2. An initialised but not build (`Dense`) keras layer (like a wrapper)
            3. A function that returns an initialised keras layer.
        kwargs: For initialisation 1. you can pass `Dense` layer kwargs

    # Input shape
        List of Atom and edge tensors of shape:
        `[(samples, max_atoms, atom_features), (samples, max_atoms, max_degrees,
          bond_features), (samples, max_atoms, max_degrees)]`
        where degrees referes to number of neighbours

    # Output shape
        New atom featuers of shape
        `(samples, max_atoms, conv_width)`

    # References
        - [Convolutional Networks on Graphs for Learning Molecular Fingerprints](https://arxiv.org/abs/1509.09292)

    '''

    def __init__(self, inner_layer_arg, activ, bias , init, **kwargs):
        # Initialise based on one of the three initialisation methods

        # Case 1: Check if inner_layer_arg is conv_width
        if isinstance(inner_layer_arg, int):
            self.conv_width = inner_layer_arg
            self.create_inner_layer_fn = lambda: layers.Dense(self.conv_width , activation = activ ,
                                                             use_bias = bias , kernel_initializer=init) ### add inputs to dense layer

        # Case 2: Check if an initialised keras layer is given
        elif isinstance(inner_layer_arg, layers.Layer):
            assert inner_layer_arg.built == False, 'When initialising with a keras layer, it cannot be built.'
            _, self.conv_width = inner_layer_arg.get_output_shape_for((None, None))
            # layer_from_config will mutate the config dict, therefore create a get fn
            self.create_inner_layer_fn = lambda: layer_from_config(dict(
                                                    class_name=inner_layer_arg.__class__.__name__,
                                                    config=inner_layer_arg.get_config()))

        # Case 3: Check if a function is provided that returns a initialised keras layer
        elif callable(inner_layer_arg):
            example_instance = inner_layer_arg()
            assert isinstance(example_instance, layers.Layer), 'When initialising with a function, the function has to return a keras layer'
            assert example_instance.built == False, 'When initialising with a keras layer, it cannot be built.'
            _, self.conv_width = example_instance.get_output_shape_for((None, None))
            self.create_inner_layer_fn = inner_layer_arg

        else:
            raise ValueError('NeuralGraphHidden has to be initialised with 1). int conv_widht, 2). a keras layer instance, or 3). a function returning a keras layer instance.')

        super(NeuralGraphHidden, self).__init__(**kwargs)

    def build(self, inputs_shape):

        # Import dimensions
        (max_atoms, max_degree, num_atom_features, num_bond_features,
         num_samples) = mol_shapes_to_dims(mol_shapes=inputs_shape)

        self.max_degree = max_degree

        # Add the dense layers (that contain trainable params)
        #   (for each degree we convolve with a different weight matrix)
        self.trainable_weights = []
        self.inner_3D_layers = []
        for degree in range(max_degree):

            # Initialise inner layer, and rename it
            inner_layer = self.create_inner_layer_fn()
            inner_layer_type = inner_layer.__class__.__name__.lower()
            inner_layer.name = self.name + '_inner_' + inner_layer_type + '_' + str(degree)

            # Initialise TimeDistributed layer wrapper in order to parallelise
            #   dense layer across atoms (3D)
            inner_3D_layer_name = self.name + '_inner_timedistributed_' + str(degree)
            inner_3D_layer = layers.TimeDistributed(inner_layer, name=inner_3D_layer_name)

            # Build the TimeDistributed layer (which will build the Dense layer)
            inner_3D_layer.build((None, max_atoms, num_atom_features+num_bond_features))

            # Store inner_3D_layer and it's weights
            self.inner_3D_layers.append(inner_3D_layer)
            self.trainable_weights += inner_3D_layer.trainable_weights

    def call(self, inputs, mask=None):
        atoms, bonds, edges = inputs

        # Import dimensions
        num_samples = atoms._keras_shape[0]
        max_atoms = atoms._keras_shape[1]
        num_atom_features = atoms._keras_shape[-1]
        num_bond_features = bonds._keras_shape[-1]

        # Create a matrix that stores for each atom, the degree it is
        atom_degrees = K.sum(tf.keras.backend.cast(K.not_equal(edges, -1),dtype = 'float32'), axis=-1, keepdims=True)

        # For each atom, look up the features of it's neighbour
        neighbour_atom_features = neighbour_lookup(atoms, edges, include_self=True)

        # Sum along degree axis to get summed neighbour features
        summed_atom_features = K.sum(neighbour_atom_features, axis=-2)

        # Sum the edge features for each atom
        summed_bond_features = K.sum(bonds, axis=-2)

        # Concatenate the summed atom and bond features
        summed_features = K.concatenate([summed_atom_features, summed_bond_features], axis=-1)

        # For each degree we convolve with a different weight matrix
        new_features_by_degree = []
        for degree in range(self.max_degree):

            # Create mask for this degree
            atom_masks_this_degree = K.cast(K.equal(atom_degrees, degree), K.floatx())

            # Multiply with hidden merge layer
            #   (use time Distributed because we are dealing with 2D input/3D for batches)
            # Add keras shape to let keras now the dimensions
            summed_features._keras_shape = (None, max_atoms, num_atom_features+num_bond_features)
            new_unmasked_features = self.inner_3D_layers[degree](summed_features)

            # Do explicit masking because TimeDistributed does not support masking
            new_masked_features = new_unmasked_features * atom_masks_this_degree

            new_features_by_degree.append(new_masked_features)

        # Finally sum the features of all atoms
        new_features = keras.layers.add(new_features_by_degree)

        return new_features

    def compute_output_shape(self, inputs_shape):

        # Import dimensions
        (max_atoms, max_degree, num_atom_features, num_bond_features,
         num_samples) = mol_shapes_to_dims(mol_shapes=inputs_shape)

        return (num_samples, max_atoms, self.conv_width)

    @classmethod
    def from_config(cls, config):
        # Use layer build function to initialise new NeuralHiddenLayer
        inner_layer_config = config.pop('inner_layer_config')
        # create_inner_layer_fn = lambda: layer_from_config(inner_layer_config.copy())
        def create_inner_layer_fn():
            return layer_from_config(deepcopy(inner_layer_config))

        layer = cls(create_inner_layer_fn, **config)
        return layer

    def get_config(self):
        config = super(NeuralGraphHidden, self).get_config()

        # Store config of (a) inner layer of the 3D wrapper
        inner_layer = self.inner_3D_layers[0].layer
        config['inner_layer_config'] = dict(config=inner_layer.get_config(),
                                            class_name=inner_layer.__class__.__name__)
        return config
		
		
		
class NeuralGraphOutput(layers.Layer):
    ''' Output Convolutional layer in a Neural Graph (as in Duvenaud et. al.,
    2015). This layer takes a graph as an input. The graph is represented as by
    three tensors.

    - The atoms tensor represents the features of the nodes.
    - The bonds tensor represents the features of the edges.
    - The edges tensor represents the connectivity (which atoms are connected to
        which)

    It returns the fingerprint vector for each sample for the given layer.

    According to the original paper, the fingerprint outputs of each hidden layer
    need to be summed in the end to come up with the final fingerprint.

    # Example
        Define the input:
        ```python
            atoms0 = Input(name='atom_inputs', shape=(max_atoms, num_atom_features))
            bonds = Input(name='bond_inputs', shape=(max_atoms, max_degree, num_bond_features))
            edges = Input(name='edge_inputs', shape=(max_atoms, max_degree), dtype='int32')
        ```

        The `NeuralGraphOutput` can be initialised in three ways:
        1. Using an integer `fp_length` and possible kwags (`Dense` layer is used)
            ```python
            fp_out = NeuralGraphOutput(fp_length, activation='relu', bias=False)([atoms0, bonds, edges])
            ```
        2. Using an initialised `Dense` layer
            ```python
            fp_out = NeuralGraphOutput(Dense(fp_length, activation='relu', bias=False))([atoms0, bonds, edges])
            ```
        3. Using a function that returns an initialised `Dense` layer
            ```python
            fp_out = NeuralGraphOutput(lambda: Dense(fp_length, activation='relu', bias=False))([atoms0, bonds, edges])
            ```

        Predict for regression:
        ```python
        main_prediction = Dense(1, activation='linear', name='main_prediction')(fp_out)
        ```

    # Arguments
        inner_layer_arg: Either:
            1. an int defining the `fp_length`, with optional kwargs for the
                inner Dense layer
            2. An initialised but not build (`Dense`) keras layer (like a wrapper)
            3. A function that returns an initialised keras layer.
        kwargs: For initialisation 1. you can pass `Dense` layer kwargs

    # Input shape
        List of Atom and edge tensors of shape:
        `[(samples, max_atoms, atom_features), (samples, max_atoms, max_degrees,
          bond_features), (samples, max_atoms, max_degrees)]`
        where degrees referes to number of neighbours

    # Output shape
        Fingerprints matrix
        `(samples, fp_length)`

    # References
        - [Convolutional Networks on Graphs for Learning Molecular Fingerprints](https://arxiv.org/abs/1509.09292)

    '''

    def __init__(self, inner_layer_arg, activ, bias , init, **kwargs):
        # Initialise based on one of the three initialisation methods

        # Case 1: Check if inner_layer_arg is fp_length
        if isinstance(inner_layer_arg, int):
            self.fp_length = inner_layer_arg
            self.create_inner_layer_fn = lambda: layers.Dense(units = self.fp_length, activation = activ ,
                                                             use_bias = bias , kernel_initializer=init) ### add inputs to dense layer

        # Case 2: Check if an initialised keras layer is given
        elif isinstance(inner_layer_arg, layers.Layer):
            assert inner_layer_arg.built == False, 'When initialising with a keras layer, it cannot be built.'
            _, self.fp_length = inner_layer_arg.get_output_shape_for((None, None))
            self.create_inner_layer_fn = lambda: inner_layer_arg

        # Case 3: Check if a function is provided that returns a initialised keras layer
        elif callable(inner_layer_arg):
            example_instance = inner_layer_arg()
            assert isinstance(example_instance, layers.Layer), 'When initialising with a function, the function has to return a keras layer'
            assert example_instance.built == False, 'When initialising with a keras layer, it cannot be built.'
            _, self.fp_length = example_instance.get_output_shape_for((None, None))
            self.create_inner_layer_fn = inner_layer_arg

        else:
            raise ValueError('NeuralGraphHidden has to be initialised with 1). int conv_widht, 2). a keras layer instance, or 3). a function returning a keras layer instance.')

        super(NeuralGraphOutput, self).__init__(**kwargs)

    def build(self, inputs_shape):

        # Import dimensions
        (max_atoms, max_degree, num_atom_features, num_bond_features,
         num_samples) = mol_shapes_to_dims(mol_shapes=inputs_shape)

        # Add the dense layer that contains the trainable parameters
        # Initialise dense layer with specified params (kwargs) and name
        inner_layer = self.create_inner_layer_fn()
        inner_layer_type = inner_layer.__class__.__name__.lower()
        inner_layer.name = self.name + '_inner_'+ inner_layer_type

        # Initialise TimeDistributed layer wrapper in order to parallelise
        #   dense layer across atoms
        inner_3D_layer_name = self.name + '_inner_timedistributed'
        self.inner_3D_layer = layers.TimeDistributed(inner_layer, name=inner_3D_layer_name)

        # Build the TimeDistributed layer (which will build the Dense layer)
        self.inner_3D_layer.build((None, max_atoms, num_atom_features+num_bond_features))

        # Store dense_3D_layer and it's weights
        self.trainable_weights = self.inner_3D_layer.trainable_weights


    def call(self, inputs, mask=None):
        
        atoms, bonds, edges = inputs

        # Import dimensions
        num_samples = atoms._keras_shape[0]
        max_atoms = atoms._keras_shape[1]
        num_atom_features = atoms._keras_shape[-1]
        num_bond_features = bonds._keras_shape[-1]

        # Create a matrix that stores for each atom, the degree it is, use it
        #   to create a general atom mask (unused atoms are 0 padded)
        # We have to use the edge vector for this, because in theory, a convolution
        #   could lead to a zero vector for an atom that is present in the molecule
        atom_degrees = K.sum(tf.keras.backend.cast(K.not_equal(edges, -1),dtype = 'float32'), axis=-1, keepdims=True)
        general_atom_mask = K.cast(K.not_equal(atom_degrees, 0), K.floatx())

        # Sum the edge features for each atom
        summed_bond_features = K.sum(bonds, axis=-2)

        # Concatenate the summed atom and bond features
        atoms_bonds_features = keras.layers.Concatenate(axis=-1)([atoms, summed_bond_features])

        # Compute fingerprint
        
        fingerprint_out_unmasked = self.inner_3D_layer(atoms_bonds_features)

        # Do explicit masking because TimeDistributed does not support masking
        fingerprint_out_masked = fingerprint_out_unmasked * general_atom_mask
        print(K.int_shape(fingerprint_out_masked))

        # Sum across all atoms
        final_fp_out = K.sum(fingerprint_out_masked, axis=-2, keepdims = False)
        print(K.int_shape(final_fp_out))

        return final_fp_out

    def compute_output_shape(self, inputs_shape):

        # Import dimensions
        (max_atoms, max_degree, num_atom_features, num_bond_features,
         num_samples) = mol_shapes_to_dims(mol_shapes=inputs_shape)

        return (num_samples, self.fp_length)

    @classmethod
    def from_config(cls, config):
        # Use layer build function to initialise new NeuralGraphOutput
        inner_layer_config = config.pop('inner_layer_config')
        create_inner_layer_fn = lambda: layer_from_config(deepcopy(inner_layer_config))

        layer = cls(create_inner_layer_fn, **config)
        return layer

    def get_config(self):
        config = super(NeuralGraphOutput, self).get_config()

        # Store config of inner layer of the 3D wrapper
        inner_layer = self.inner_3D_layer.layer
        config['inner_layer_config'] = dict(config=inner_layer.get_config(),
                                            class_name=inner_layer.__class__.__name__)
        return config