import json
import pytest
import numpy as np

from pssm import MSA
from pathlib import Path


@pytest.fixture(scope="module")
def json_data():
    test_json = 'pssm_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


def get_sequences():
    test_json = 'pssm_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data['msa_sequences']


def get_bg():
    test_json = 'pssm_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data['bg_matrix']


@pytest.fixture(scope="module")
def msa_sequences(json_data):
    return json_data['msa_sequences']


@pytest.fixture(scope="module")
def invalid_msa(json_data):
    return json_data['invalid_msa']


@pytest.fixture(scope="module")
def primary_seq(json_data):
    return json_data['primary_seq']


@pytest.fixture(scope="module")
def msa_size(json_data):
    return tuple(json_data['msa_size'])


@pytest.fixture(scope="module")
def sequence_weights(json_data):
    return np.array(json_data['sequence_weights'])


@pytest.fixture(scope="module")
def num_observations(json_data):
    return np.float64(json_data['num_observations'])


@pytest.fixture(scope="module")
def bg_matrix(json_data):
    return np.array(json_data['bg_matrix'])


@pytest.fixture(scope="module")
def pssm_01(json_data):
    return np.array(json_data['pssm_01'])


@pytest.fixture(scope="module")
def pssm_02(json_data):
    return np.array(json_data['pssm_02'])


@pytest.fixture(scope="module")
def pssm_03(json_data):
    return np.array(json_data['pssm_03'])


@pytest.fixture(scope="module")
def pssm_04(json_data):
    return np.array(json_data['pssm_04'])


@pytest.fixture(scope="module")
def pssm_05(json_data):
    return np.array(json_data['pssm_05'])


@pytest.fixture(scope="module")
def pssm_06(json_data):
    return np.array(json_data['pssm_06'])


@pytest.fixture(scope="module")
def pssm_07(json_data):
    return np.array(json_data['pssm_07'])


@pytest.fixture(scope="module")
def pssm_08(json_data):
    return np.array(json_data['pssm_08'])


@pytest.fixture(scope="module")
def pssm_09(json_data):
    return np.array(json_data['pssm_09'])


@pytest.fixture(scope="module")
def pssm_10(json_data):
    return np.array(json_data['pssm_10'])


msa = None


def check_init():
    assert msa is not None, 'MSA initialization failed.'


def check_error(sequences):
    try:
        msa = MSA(sequences)
    except TypeError:
        return True
    except Exception:
        return False
    return False


def test_pssm_raise_error(msa_sequences, invalid_msa):
    global msa

    try:
        msa = MSA(msa_sequences)
    except Exception:
        msa = None

    check_init()

    passed = all([check_error(inv_msa) for inv_msa in invalid_msa])

    assert passed, 'Failed to raise TypeError for invalid MSA.'


def test_pssm_get_size(msa_size):
    check_init()

    passed = True

    try:
        size = msa.get_size()
    except Exception:
        passed = False

    assert passed, 'Error while testing get_size().'

    passed = (size == msa_size)

    assert passed, 'Incorrect MSA size.'


def test_pssm_get_primary_sequence(primary_seq):
    check_init()

    passed = True

    try:
        seq = msa.get_primary_sequence()
    except Exception:
        passed = False

    assert passed, 'Error while testing get_primary_sequence().'

    passed = (seq == primary_seq)

    assert passed, 'Incorrect primary sequence.'


def test_pssm_get_sequence_weights(sequence_weights):
    check_init()

    passed = True

    try:
        weights = msa.get_sequence_weights()
    except Exception:
        passed = False

    assert passed, 'Error while testing get_sequence_weights().'

    try:
        passed = all(np.isclose(weights, sequence_weights, atol=1e-8, rtol=0))
    except Exception:
        passed = False

    assert passed, 'Incorrect sequence weights.'


def test_pssm_get_number_of_observations(num_observations):
    check_init()

    passed = True

    try:
        num_obs = msa.get_number_of_observations()
    except Exception:
        passed = False

    assert passed, 'Error while testing get_number_of_observations().'

    try:
        passed = np.isclose(num_obs, num_observations, atol=1e-8, rtol=0)
    except Exception:
        passed = False

    assert passed, 'Incorrect number of independent observations.'


def check_pssm(pssm):
    passed = ((type(pssm) == np.ndarray) and (pssm.dtype == np.int64))

    assert passed, 'Return value not a numpy.ndarray of dtype numpy.int64.'


def test_pssm_get_pssm_basic(pssm_01):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm()
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm().'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_01)

    assert passed, 'Incorrect PSSM (basic).'


def test_pssm_get_pssm_with_bg_matrix(pssm_02, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm(bg_matrix).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_02)

    assert passed, 'Incorrect PSSM (bg_matrix).'


def test_pssm_get_pssm_with_redistribute_gaps(pssm_03):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(redistribute_gaps=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm(redistribute_gaps).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_03)

    assert passed, 'Incorrect PSSM (redistribute_gaps).'


def test_pssm_get_pssm_with_sequence_weights(pssm_04):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(use_sequence_weights=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm(use_sequence_weights).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_04)

    assert passed, 'Incorrect PSSM (use_sequence_weights).'


def test_pssm_get_pssm_with_bg_matrix_and_redistribute_gaps(pssm_05, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix, redistribute_gaps=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm' \
                   '(bg_matrix, redistribute_gaps).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_05)

    assert passed, 'Incorrect PSSM (bg_matrix, redistribute_gaps).'


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights(pssm_06, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix, use_sequence_weights=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm' \
                   '(bg_matrix, use_sequence_weights).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_06)

    assert passed, 'Incorrect PSSM (bg_matrix, use_sequence_weights).'


def test_pssm_get_pssm_with_pseudocounts(pssm_07, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(add_pseudocounts=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm(add_pseudocounts).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_07)

    assert passed, 'Incorrect PSSM (add_pseudocounts).'


def test_pssm_get_pssm_with_bg_matrix_and_redistribute_gaps_and_pseudocounts(pssm_08, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            redistribute_gaps=True,
                            add_pseudocounts=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm' \
                   '(bg_matrix, redistribute_gaps, add_pseudocounts).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_08)

    assert passed, 'Incorrect PSSM ' \
                   '(bg_matrix, redistribute_gaps, add_pseudocounts).'


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_pseudocounts(pssm_09, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            use_sequence_weights=True,
                            add_pseudocounts=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm' \
                   '(bg_matrix, use_sequence_weights, add_pseudocounts).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_09)

    assert passed, 'Incorrect PSSM ' \
                   '(bg_matrix, use_sequence_weights, add_pseudocounts).'


def test_pssm_get_pssm_with_bg_matrix_and_sequence_weights_and_redistribute_gaps_and_pseudocounts(pssm_10, bg_matrix):
    check_init()

    passed = True

    try:
        pssm = msa.get_pssm(bg_matrix=bg_matrix,
                            redistribute_gaps=True,
                            use_sequence_weights=True,
                            add_pseudocounts=True)
    except Exception:
        passed = False

    assert passed, 'Error while testing get_pssm' \
                   '(bg_matrix, redistribute_gaps, ' \
                   'use_sequence_weights, add_pseudocounts).'

    check_pssm(pssm)

    passed = np.array_equal(pssm, pssm_10)

    assert passed, 'Incorrect PSSM ' \
                   '(bg_matrix, redistribute_gaps, ' \
                   'use_sequence_weights, add_pseudocounts).'
