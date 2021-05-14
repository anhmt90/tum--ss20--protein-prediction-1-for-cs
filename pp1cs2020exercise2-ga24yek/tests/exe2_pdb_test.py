import json
import pytest
import numpy as np

from pathlib import Path
from exe2_pdb import PDB_Parser


@pytest.fixture(scope="module")
def relative_path():
    return Path(__file__).parent


@pytest.fixture(scope="module")
def json_data(relative_path):
    test_json = 'exe2_pdb_test.json'

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def filename(relative_path, json_data):
    return Path(relative_path, json_data['filename'])


@pytest.fixture(scope="module")
def bfactors_1(relative_path, json_data):
    return np.load(Path(relative_path, json_data["bfactors_filename_1"]))


@pytest.fixture(scope="module")
def bfactors_2(relative_path, json_data):
    return np.load(Path(relative_path, json_data["bfactors_filename_2"]))


@pytest.fixture(scope="module")
def bfactors_3(relative_path, json_data):
    return np.load(Path(relative_path, json_data["bfactors_filename_3"]))


@pytest.fixture(scope="module")
def bfactors_4(relative_path, json_data):
    return np.load(Path(relative_path, json_data["bfactors_filename_4"]))


@pytest.fixture(scope="module")
def bfactors_5(relative_path, json_data):
    return np.load(Path(relative_path, json_data["bfactors_filename_5"]))


@pytest.fixture(scope="module")
def contact_map_1(relative_path, json_data):
    return np.load(Path(relative_path, json_data["contact_map_filename_1"]))


@pytest.fixture(scope="module")
def contact_map_2(relative_path, json_data):
    return np.load(Path(relative_path, json_data["contact_map_filename_2"]))


@pytest.fixture(scope="module")
def contact_map_3(relative_path, json_data):
    return np.load(Path(relative_path, json_data["contact_map_filename_3"]))


@pytest.fixture(scope="module")
def contact_map_4(relative_path, json_data):
    return np.load(Path(relative_path, json_data["contact_map_filename_4"]))


@pytest.fixture(scope="module")
def contact_map_5(relative_path, json_data):
    return np.load(Path(relative_path, json_data["contact_map_filename_5"]))


pdb_parser = None


def check_init():
    assert pdb_parser is not None, 'Failed to initiliaze PDB_Parser.'


def test_pdb_get_number_of_chains(json_data, filename):
    global pdb_parser

    passed = True

    try:
        pdb_parser = PDB_Parser(filename)
    except Exception:
        passed = False
        pdb_parser = None

    assert passed, 'Error while creating PDB_Parser.'

    check_init()

    try:
        num_chains = pdb_parser.get_number_of_chains()
    except Exception:
        passed = False

    hint = 'General structure of Bio.PDB: Structure -> Model -> Chains -> Residues -> Atoms.'

    assert passed, 'Error in PDB_Parser.get_number_of_chains().'

    passed = (json_data['number_of_chains'] == num_chains)

    assert passed, f'Incorrect number of chains. {hint}'


def test_pdb_get_sequence(json_data):
    check_init()

    passed = True

    try:
        seq = pdb_parser.get_sequence(json_data['sequence_chain'])
    except Exception:
        passed = False

    assert passed, 'Error in PDB_Parser.get_sequence().'

    passed = (json_data['sequence'] == seq)

    assert passed, 'Incorrect sequence.'


def test_pdb_get_number_of_water_molecules(json_data):
    check_init()

    passed = True

    try:
        chain = json_data['water_molecules_chain']
        num_water = pdb_parser.get_number_of_water_molecules(chain)
    except Exception:
        passed = False

    hint = 'Water molecules are called HOH and are part of Bio.PDBs chains.'

    assert passed, 'Error in PDB_Parser.get_number_of_water_molecules().'

    passed = (json_data['number_of_water_molecules'] == num_water)

    assert passed, f'Incorrect number of water molecules. {hint}'


def test_pdb_get_ca_distance_same_chain(json_data):
    check_init()

    passed = True

    try:
        ca_dis = pdb_parser.get_ca_distance(
            json_data["ca_distance_same_chains"]["chain_1"],
            json_data["ca_distance_same_chains"]["id_1"],
            json_data["ca_distance_same_chains"]["chain_2"],
            json_data["ca_distance_same_chains"]["id_2"]
        )
    except Exception:
        passed = False

    assert passed, 'Error in PDB_Parser.get_ca_distance().'

    passed = (json_data["ca_distance_same_chains"]["result"] == ca_dis)

    assert passed, 'Incorrect C-alpha distance (same chain).'


def test_pdb_get_ca_distance_different_chains(json_data):
    check_init()

    passed = True

    try:
        ca_dis = pdb_parser.get_ca_distance(
            json_data["ca_distance_diff_chains"]["chain_1"],
            json_data["ca_distance_diff_chains"]["id_1"],
            json_data["ca_distance_diff_chains"]["chain_2"],
            json_data["ca_distance_diff_chains"]["id_2"]
        )
    except Exception:
        passed = False

    assert passed, 'Error in PDB_Parser.get_ca_distance().'

    passed = (json_data["ca_distance_diff_chains"]["result"] == ca_dis)

    assert passed, 'Incorrect C-alpha distance (different chains).'


def check_bfactors(chain, bfactors):
    check_init()

    passed = True

    try:
        bf = pdb_parser.get_bfactors(chain)
    except Exception:
        passed = False

    assert passed, 'Error in PDB_Parser.get_bfactors().'

    passed = (type(bf) == np.ndarray and bf.dtype == np.int64)

    assert passed, 'B-Factors are not a numpy.ndarray of dtype numpy.int64.'

    passed = (bfactors.shape == bf.shape)

    assert passed, 'The shape of the numpy array is incorrect.'

    passed = np.all((bfactors == bf) | (np.isnan(bfactors) & np.isnan(bf)))

    assert passed, 'Incorrect B-Factors.'


def test_pdb_get_bfactors_1(json_data, bfactors_1):
    check_bfactors(json_data['bfactors_chain_1'], bfactors_1)


def test_pdb_get_bfactors_2(json_data, bfactors_2):
    check_bfactors(json_data['bfactors_chain_2'], bfactors_2)


def test_pdb_get_bfactors_3(json_data, bfactors_3):
    check_bfactors(json_data['bfactors_chain_3'], bfactors_3)


def test_pdb_get_bfactors_4(json_data, bfactors_4):
    check_bfactors(json_data['bfactors_chain_4'], bfactors_4)


def test_pdb_get_bfactors_5(json_data, bfactors_5):
    check_bfactors(json_data['bfactors_chain_5'], bfactors_5)


def check_contact_map(chain, contact_map):
    check_init()

    passed = True

    try:
        c_map = pdb_parser.get_contact_map(chain)
    except Exception:
        passed = False

    assert passed, 'Error in PDB_Parser.get_contact_map().'

    passed = (type(c_map) == np.ndarray and c_map.dtype == np.int64)

    assert passed, 'Contact map is not a numpy.ndarray of dtype numpy.int64.'

    passed = (contact_map.shape == c_map.shape)

    assert passed, f'The shape of the numpy array is incorrect.'

    passed = np.all((contact_map == c_map)
                    | (np.isnan(contact_map) & np.isnan(c_map)))

    assert passed, 'Incorrect contact map.'


def test_pdb_get_contact_map_1(json_data, contact_map_1):
    check_contact_map(json_data['contact_map_chain_1'], contact_map_1)


def test_pdb_get_contact_map_2(json_data, contact_map_2):
    check_contact_map(json_data['contact_map_chain_2'], contact_map_2)


def test_pdb_get_contact_map_3(json_data, contact_map_3):
    check_contact_map(json_data['contact_map_chain_3'], contact_map_3)


def test_pdb_get_contact_map_4(json_data, contact_map_4):
    check_contact_map(json_data['contact_map_chain_4'], contact_map_4)


def test_pdb_get_contact_map_5(json_data, contact_map_5):
    check_contact_map(json_data['contact_map_chain_5'], contact_map_5)
