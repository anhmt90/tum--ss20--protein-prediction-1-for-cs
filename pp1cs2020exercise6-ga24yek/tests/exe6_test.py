# -*- coding: utf-8 -*-
"""
Created on Thu Jun  7 16:11:18 2018

@author: Michael
"""

import pytest
import json

import os
import math
import numpy as np

from exe6_perceptron import Perceptron


############ HELPER FUNCTIONS ##################
@pytest.fixture(scope="module")
def relative_path():
    return os.path.dirname(__file__)


@pytest.fixture(scope="module")
def json_data(relative_path):
    with open(os.path.join(relative_path, 'exe6_test.json')) as json_file:
        json_data = json.load(json_file)
    return json_data


@pytest.fixture(scope="module")
def bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["bias"]))


@pytest.fixture(scope="module")
def hinge_loss(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["hinge"]))


@pytest.fixture(scope="module")
def delta_hinge(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["delta_hinge"]))


@pytest.fixture(scope="module")
def l2_loss(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["l2_loss"]))


@pytest.fixture(scope="module")
def delta_l2(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["delta_l2"]))


@pytest.fixture(scope="module")
def sigmoid(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["sigmoid"]))


@pytest.fixture(scope="module")
def perceptron(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["perceptron"]))

@pytest.fixture(scope="module")
def perceptron_bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["perceptron_bias"]))


@pytest.fixture(scope="module")
def multiperceptron_bias(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["multiperceptron_bias"]))


@pytest.fixture(scope="module")
def multiperceptron_bias_nonlin(relative_path, json_data):
    return np.load(os.path.join(relative_path, json_data["multiperceptron_bias_nonlin"]))

############ INIT STUDENT PERCEPTRON ######################
@pytest.fixture(scope="module")
def student_perceptron(relative_path, json_data):
    # learning rate, number of epochs and random seed for single perceptrons
    LEARNING_RATE = json_data['parameters']['learning_rate']
    NEPOCHS       = json_data['parameters']['nepochs']
    SEED          = json_data['parameters']['seed']

    try:
        perc = Perceptron( LEARNING_RATE, NEPOCHS, SEED )
    except:
        assert False, 'Something went wrong with initializing your Perceptron!'
    return perc
    

############ TESTS ##################
def test_adding_bias_term( bias, student_perceptron,  ):
    passed = True 

    # create 1D test array and 2D test array
    test_array_1D   = bias['bias_1d_in']
    test_array_2D   = bias['bias_2d_in']
    
    try:
        student_answer_1D  = student_perceptron._add_bias( test_array_1D )
        student_answer_2D  = student_perceptron._add_bias( test_array_2D )
    except:
        passed = False

    assert passed, 'Error in bias adding procedure.'
    assert student_answer_1D is not None, "No array was returned by the _add_bias function."
    assert student_answer_2D is not None, "No array was returned by the _add_bias function."
    
    try:
        correct_answer_1D = np.all( np.isclose( bias['bias_1d_out'], student_answer_1D ) )
        correct_answer_2D = np.all( np.isclose( bias['bias_2d_out'], student_answer_2D ) )
    except:
        passed = False

    assert passed, 'Something with your bias return value went wrong.'
    assert correct_answer_1D, ( "When adding a bias term to a 1D array you have to " + 
                                "add a 1 to the end of the array." )
    
    assert correct_answer_2D, ( "When adding a bias term to a 2D array you have to " +
                                "add a 1 to the end of each row in the array.")
    
    
def test_hinge_loss( hinge_loss, student_perceptron,  ):
    passed = True

    for _, data in hinge_loss.items():
        y      = data[0]
        y_pred = data[1]

        try:
            student_answer_hinge_loss  = student_perceptron._hinge_loss( y, y_pred )
        except:
            passed = False
        
        assert passed, 'Your hinge loss function produces an error.'
        assert student_answer_hinge_loss is not None, ( 
                            "You returned None instead of the hinge loss." )
        try:
            correct_answer = math.isclose( data[2], student_answer_hinge_loss)
        except:
            passed = False
        
        assert passed, 'Your hinge loss does not seem to be well formated.' 
        assert correct_answer, ( "Your definition of hinge loss is not correct. " + 
                             "Please keep in mind that this is the loss and not the derivative! " +
                             "Also keep in mind that groundtruth labels and predicted labels are " +
                             "within [-1, 1], not within [0, 1].")

    
def test_delta_hinge( delta_hinge, student_perceptron,  ):
    passed = True

    for _, data in delta_hinge.items():
        y      = data[0]
        y_pred = data[1]

        try:
            student_answer_delta_hinge  = student_perceptron._delta_hinge( y, y_pred )
        except:
            passed = False
        
        assert passed, 'Your delta hinge function seems to produce an error.'
        assert student_answer_delta_hinge is not None, (
                    "You returned None instead of the derivative hinge loss." )

        try:
            correct_answer = math.isclose( data[2], student_answer_delta_hinge)
        except:
            passed = False

        assert passed, 'Your delta hinge does not seem to be well formated.' 
        assert correct_answer, ( "Your definition of the derivative of the hinge loss is not correct.")

    
def test_l2_loss( l2_loss, student_perceptron , ):
    passed = True 

    for _, data in l2_loss.items():
        y      = data[0]
        y_pred = data[1]

        try: 
            student_answer_l2_loss = student_perceptron._l2_loss( y, y_pred )
        except:
            passed = False
        
        assert passed, 'Your l2 loss seems to produce an error.'
        assert student_answer_l2_loss is not None, ( 
                        "You returned None instead of the l2 loss." )
        try:
            correct_answer = math.isclose( data[2], student_answer_l2_loss)
        except:
            passed = False

        assert passed, 'Your l2 loss does not seem to be well formated.' 
        assert correct_answer, ( "Your definition of l2 loss is not correct. " + 
                             "Please keep in mind that this is the loss and not the derivative.")


def test_delta_l2( delta_l2, student_perceptron,  ):
    passed = True

    for _, data in delta_l2.items():
        y      = data[0]
        y_pred = data[1]
        try:
            student_answer_delta_l2 = student_perceptron._delta_l2( y, y_pred )
        except:
            passed = False

        assert passed, 'Your delta l2 function seems to produce an error.'
        assert student_answer_delta_l2 is not None, (
                "You returned None instead of the correct derivative of the l2 loss." )

        try:
            correct_answer = math.isclose( data[2], student_answer_delta_l2 )
        except:
            passed = False
        assert passed, 'Your delta l2 does not seem to be well formated.' 
        assert correct_answer, ( "Your definition of the derivative of the l2 loss is not correct.")
        
    
def test_sigmoid( sigmoid, student_perceptron,  ):
    passed = True

    try:
        student_answer_sigmoid = student_perceptron._sigmoid( sigmoid[0,:] )
    except:
        passed = False

    assert passed, 'Your sigmoid function seems to produce an error.'   
    assert student_answer_sigmoid is not None,  (
                    "You returned None instead of the correct sigmoid." )
    try:
        correct_answer = np.all( np.isclose( sigmoid[1,:], student_answer_sigmoid) )
    except:
        passed = False
    assert passed, 'Your sigmoid does not seem to be well formated.' 
    assert correct_answer, ( "Your definition of sigmoid is not correct.")
    
    
############################################################################## 
################################ 2 POINTS ####################################
    
def test_single_perceptron( perceptron, student_perceptron,  ):
    passed = True

    try:
        student  = student_perceptron.single_perceptron()
    except:
        passed = False
    
    assert passed, 'Your delta single_perceptron function seems to produce an error.'
    assert student is not None, ( 
                    "You did not return any weights for the single_perceptron." )
    
    try:
        correct_answer = np.all( np.isclose( perceptron, student ) )
    except:
        passed = False
    
    assert passed, 'Your single_perceptron does not seem to be well formated.' 
    assert correct_answer, ( "Your single_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs."  +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: hinge loss) and the correct targets (OR gate)." )

def test_single_perceptron1( perceptron, student_perceptron,  ):
    passed = True

    try:
        student  = student_perceptron.single_perceptron()
    except:
        passed = False
    
    assert passed, 'Your single_perceptron function seems to produce an error.'
    assert student is not None, ( 
                    "You did not return any weights for the single_perceptron." )
    
    try:
        correct_answer = np.all( np.isclose( perceptron, student ) )
    except:
        passed = False
    
    assert passed, 'Your single_perceptron does not seem to be well formated.' 
    assert correct_answer, ( "Your single_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs."  +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: hinge loss) and the correct targets (OR gate)." )
    
############################################################################## 
################################ 2 POINTS ####################################
def test_single_perceptron_with_bias( perceptron_bias, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.single_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your single_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
                "You did not return any weights for the single_perceptron with a bias term." +
                "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( perceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your single_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your single_perceptron_with_bias did not return the " +
                                "correct weights based on the given learning rate " +
                                "and number of epochs. Double check whether you "   +
                                "are using the correct loss function ( here: hinge loss) " + 
                                "and the correct targets (OR gate)." )


def test_single_perceptron_with_bias1( perceptron_bias, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.single_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your single_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
                "You did not return any weights for the single_perceptron with a bias term." +
                "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( perceptron_bias, student ) )
    except:
        passed = False
        
    assert passed, 'Your single_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your single_perceptron_with_bias did not return the " +
                                "correct weights based on the given learning rate " +
                                "and number of epochs. Double check whether you "   +
                                "are using the correct loss function ( here: hinge loss) " + 
                                "and the correct targets (OR gate)." )
    
    
    
##################### START TESTING MULTILAYER PERCEPTRONS ###################
################################ 5 POINTS ####################################
    
def test_multi_perceptron_with_bias( multiperceptron_bias, student_perceptron,  ):
    passed = True
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )


def test_multi_perceptron_with_bias1( multiperceptron_bias, student_perceptron,  ):
    passed = True
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )
    
    
def test_multi_perceptron_with_bias2( multiperceptron_bias, student_perceptron,  ):
    passed = True
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )
    
    
    
def test_multi_perceptron_with_bias3( multiperceptron_bias, student_perceptron,  ):
    passed = True
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )



def test_multi_perceptron_with_bias4( multiperceptron_bias, student_perceptron,  ):
    passed = True
    try:
        student = student_perceptron.multi_perceptron_with_bias()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )
    
    
############################################################################## 
################################ 5 POINTS ####################################
def test_multi_perceptron_with_bias_and_nonlinearity( multiperceptron_bias_nonlin, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias_nonlin, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )

def test_multi_perceptron_with_bias_and_nonlinearity1( multiperceptron_bias_nonlin, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias_nonlin, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )

def test_multi_perceptron_with_bias_and_nonlinearity2( multiperceptron_bias_nonlin, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias_nonlin, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )

def test_multi_perceptron_with_bias_and_nonlinearity3( multiperceptron_bias_nonlin, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias_nonlin, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )

def test_multi_perceptron_with_bias_and_nonlinearity4( multiperceptron_bias_nonlin, student_perceptron,  ):
    passed = True

    try:
        student = student_perceptron.multi_perceptron_with_bias_and_nonlinearity()
    except:
        passed = False
        student = None

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity function seems to produce an error.'
    assert student is not None, ( 
            "You did not return any weights for the multi_perceptron with a bias term. "+
            "Probably, you did not implement the _add_bias function, yet." )
    try:
        correct_answer = np.all( np.isclose( multiperceptron_bias_nonlin, student ) )
    except:
        passed = False

    assert passed, 'Your multi_perceptron_with_bias_and_nonlinearity does not seem to be well formated.' 
    assert correct_answer, ( "Your multi_perceptron did not return the correct weights " +
                                "based on the given learning rate and number of epochs. " +
                                "Double check whether you are using the correct loss "    + 
                                " function ( here: l2 loss) and the correct targets (XOR gate)." )
    

    
