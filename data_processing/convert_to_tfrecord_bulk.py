#!/usr/bin/python

# imports
import sys
import re
import tensorflow as tf
import os

# Constants
NUM_DIMENSIONS = 3

# Accesory functions for dealing with TF Example and SequenceExample
_example = tf.train.Example
_sequence_example = tf.train.SequenceExample
_feature = tf.train.Feature
_features = lambda d: tf.train.Features(feature=d)
_feature_list = lambda l: tf.train.FeatureList(feature=l)
_feature_lists = lambda d: tf.train.FeatureLists(feature_list=d)
_bytes_feature = lambda v: _feature(bytes_list=tf.train.BytesList(value=v))
_int64_feature = lambda v: _feature(int64_list=tf.train.Int64List(value=v))
_float_feature = lambda v: _feature(float_list=tf.train.FloatList(value=v))

# Functions for conversion from Text-based ProteinNet files to TFRecords
_aa_dict = {'A': '0', 'C': '1', 'D': '2', 'E': '3', 'F': '4', 'G': '5', 'H': '6', 'I': '7', 'K': '8', 'L': '9',
            'M': '10', 'N': '11', 'P': '12', 'Q': '13', 'R': '14', 'S': '15', 'T': '16', 'V': '17', 'W': '18',
            'Y': '19'}
_dssp_dict = {'L': '0', 'H': '1', 'B': '2', 'E': '3', 'G': '4', 'I': '5', 'T': '6', 'S': '7'}
_mask_dict = {'-': '0', '+': '1'}


class switch(object):
    """Switch statement for Python, based on recipe from Python Cookbook."""

    def __init__(self, value):
        self.value = value
        self.fall = False

    def __iter__(self):
        """Return the match method once, then stop"""
        yield self.match
        raise StopIteration

    def match(self, *args):
        """Indicate whether or not to enter a case suite"""
        if self.fall or not args:
            return True
        elif self.value in args:  # changed for v1.5
            self.fall = True
            return True
        else:
            return False


def letter_to_num(string, dict_):
    """ Convert string of letters to list of ints """
    patt = re.compile('[' + ''.join(dict_.keys()) + ']')
    num_string = patt.sub(lambda m: dict_[m.group(0)] + ' ', string)
    num = [int(i) for i in num_string.split()]
    return num


def read_record(file_, num_evo_entries):
    """ Read a Mathematica protein record from file and convert into dict. """

    dict_ = {}
    notDone = True #go until "\n" is found
    next_line = file_.readline()
    i = 0
    lastKey = ""
    while next_line:
        if next_line == "\n":
            print ("DONE")
        elif "[" in next_line and "]" in next_line:
            lastKey = next_line.strip().replace("[", "").replace("]","").lower() #remove end line
            print (lastKey)
        else:
            if lastKey == "id":
                dict_[lastKey] = next_line.strip()
            elif lastKey == "primary":
                dict_[lastKey] = primary = letter_to_num(next_line.strip(), _aa_dict)
            elif lastKey == "evolutionary":
                evolutionary = []
                for residue in range(num_evo_entries): evolutionary.append(
                    [float(step) for step in file_.readline().split()])
                dict_['evolutionary'] = evolutionary
            elif lastKey == "secondary":
                secondary = letter_to_num(next_line.strip(), _dssp_dict)
                dict_['secondary'] = secondary
            elif lastKey == "tertiary":
                tertiary = []
                for axis in range(NUM_DIMENSIONS): tertiary.append(
                    [float(coord) for coord in file_.readline().split()])
                dict_['tertiary'] =  tertiary
            elif lastKey == "mask":
                mask = letter_to_num(file_.readline().strip(), _mask_dict)
                dict_.update({'mask': mask})
            else:
                print ("UNRECOGNIZED!")
        next_line = file_.readline()
    print (dict_)
    return dict_


def dict_to_tfrecord(dict_):
    """ Convert protein dict into TFRecord. """
    bytes = dict_['id'].encode('utf-8')
    print(bytes)
    id_ = _bytes_feature( [bytes])

    feature_lists_dict = {}
    feature_lists_dict.update({'primary': _feature_list([_int64_feature([aa]) for aa in dict_['primary']])})

    if 'evolutionary' in dict_:
        feature_lists_dict.update(
            {'evolutionary': _feature_list([_float_feature(list(step)) for step in zip(*dict_['evolutionary'])])})

    if 'secondary' in dict_:
        feature_lists_dict.update({'secondary': _feature_list([_int64_feature([dssp]) for dssp in dict_['secondary']])})

    if 'tertiary' in dict_:
        feature_lists_dict.update(
            {'tertiary': _feature_list([_float_feature(list(coord)) for coord in zip(*dict_['tertiary'])])})

    if 'mask' in dict_:
        feature_lists_dict.update({'mask': _feature_list([_float_feature([step]) for step in dict_['mask']])})

    record = _sequence_example(context=_features({'id': id_}), feature_lists=_feature_lists(feature_lists_dict))

    return record

def processOne(input_path, output_path, num_evo_entries = 42):
    input_file = open(input_path, 'r')
    print (input_file)
    output_file = tf.compat.v1.python_io.TFRecordWriter(output_path)
    while True:
        dict_ = read_record(input_file, num_evo_entries)
        print (dict_)
        if dict_:
            tfrecord_serialized = dict_to_tfrecord(dict_).SerializeToString()
            print ("serialized")
            output_file.write(tfrecord_serialized)
            print ("wrote out")
        else:
            input_file.close()
            output_file.close()
            break
    print ("processed one!")

# main. accepts three command-line arguments: input file, output file, and the number of entries in evo profiles
if __name__ == '__main__':
    inputsDir = "../pssms/"#sys.argv[1] #directory with proteinnet records
    outputsDir = "../tfRecords/"#sys.argv[2] #directory to put tfrecords
    #for each proteinnet in the folder, process and output a tfrecord
    for file in os.listdir(inputsDir):
        if '.proteinnet' in file:
            outputFile = outputsDir + file.replace(".proteinnet","") + ".tfrecord"
            print("input: ", file, " output: ", outputFile)
            processOne(inputsDir + file, outputFile)


