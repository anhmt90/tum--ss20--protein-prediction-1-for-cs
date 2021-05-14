import json
import pytest

from pathlib import Path
from exe2_swissprot import SwissProt_Parser


@pytest.fixture(scope="module")
def relative_path():
    return Path(__file__).parent


@pytest.fixture(scope="module")
def json_data(relative_path):
    test_json = 'exe2_swissprot_test.json'

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def filename(relative_path, json_data):
    return Path(relative_path, json_data['filename'])


@pytest.fixture(scope="module")
def identifier(json_data):
    return json_data['identifier']


@pytest.fixture(scope="module")
def sequence_length(json_data):
    return json_data['sequence_length']


@pytest.fixture(scope="module")
def organism(json_data):
    return json_data['organism']


@pytest.fixture(scope="module")
def localization(json_data):
    return set(json_data['localization'])


@pytest.fixture(scope="module")
def pdb_support(json_data):
    return set(json_data['pdb_support'])


sp_parser = None


def check_init():
    assert sp_parser is not None, 'Failed to initiliaze SwissProt_Parser.'


def test_sp_get_sp_identifier(filename, identifier):
    global sp_parser

    passed = True

    try:
        sp_parser = SwissProt_Parser(filename)
    except Exception:
        passed = False
        sp_parser = None

    assert passed, 'Error while creating SwissProt_Parser.'

    check_init()

    try:
        sp_id = sp_parser.get_sp_identifier()
    except Exception:
        passed = False

    assert passed, 'Error in SwissProt_Parser.get_sp_identifier().'

    passed = (identifier == sp_id)

    assert passed, 'Incorrect SwissProt identifier.'


def test_sp_get_sp_sequence_length(sequence_length):
    check_init()

    passed = True

    try:
        seq_length = sp_parser.get_sp_sequence_length()
    except Exception:
        passed = False

    assert passed, 'Error in SwissProt_Parser.get_sp_sequence_length().'

    passed = (sequence_length == seq_length)

    assert passed, 'Incorrect sequence length.'


def test_sp_get_organism(organism):
    check_init()

    passed = True

    try:
        org = sp_parser.get_organism()
    except Exception:
        passed = False

    hint = 'Check the field [organism].'

    assert passed, 'Error in SwissProt_Parser.get_organism().'

    passed = (organism == org)

    assert passed, f'Incorrect organism. {hint}'


def test_sp_get_localization(localization):
    check_init()

    passed = True

    try:
        loc = sp_parser.get_localization()
    except Exception:
        passed = False

    hint = 'Check the field [comment_subcellularlocation_location].'

    assert passed, 'Error in SwissProt_Parser.get_localization().'

    passed = (len(localization) == len(loc))

    assert passed, f'Incorrect number of localizations. {hint}'

    passed = (localization == set(loc))

    assert passed, f'Incorrect localization(s). {hint}'


def test_sp_get_pdb_support(pdb_support):
    check_init()

    passed = True

    try:
        pdb = sp_parser.get_pdb_support()
    except Exception:
        passed = False

    hint = 'Check the field [dbxrefs].'

    assert passed, 'Error in SwissProt_Parser.get_pdb_support().'

    passed = (len(pdb_support) == len(pdb))

    assert passed, f'Incorrect number of PDB IDs. {hint}'

    passed = (pdb_support == set(pdb))

    assert passed, f'Incorrect PDB ID(s). {hint}'
