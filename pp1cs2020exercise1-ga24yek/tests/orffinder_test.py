import json
import pytest

from math import ceil
from pathlib import Path
from orffinder import get_orfs


@pytest.fixture(scope="module")
def json_data():
    test_json = 'orffinder_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def sequence(json_data):
    return json_data['sequence']


@pytest.fixture(scope="module")
def invalid_genome(json_data):
    return json_data['invalid_genome']


@pytest.fixture(scope="module")
def min_num_aa(json_data):
    return json_data['min_num_aa']


@pytest.fixture(scope="module")
def orf_list(json_data):
    orfs = [tuple(e) for e in json_data['orf_list']]
    return set(orfs)


valid_orfs = 0
total_orfs = 999999
student_orfs = None


def check_init():
    assert student_orfs is not None, 'Did not pass test_orffinder_get_orfs().'


def test_orffinder_get_orfs(sequence, min_num_aa):
    global student_orfs

    passed = True
    student_orfs = None

    try:
        orfs = get_orfs(sequence, min_num_aa)
    except Exception:
        passed = False
    finally:
        student_orfs = None

    assert passed, 'Error in get_orfs().'

    passed = isinstance(orfs, list)

    assert passed, 'Return type is not a list.'

    passed = all([isinstance(t, tuple) for t in orfs])

    assert passed, 'List does not contain only tuples.'

    passed = all([(len(t) == 4) for t in orfs])

    assert passed, 'Not all tuples contain exactly 4 elements.'

    for p1, p2, seq, flag in orfs:
        passed = all([isinstance(p1, int),
                      isinstance(p2, int),
                      isinstance(seq, str),
                      isinstance(flag, bool)])

        assert passed, 'Tuple malformed. Expected: (int, int, str, bool).'

    student_orfs = orfs


def test_orffinder_raise_error(invalid_genome, min_num_aa):
    passed = False

    check_init()

    try:
        get_orfs(invalid_genome, min_num_aa)
    except TypeError:
        passed = True
    except Exception:
        passed = False

    assert passed, 'Failed to raise TypeError for invalid genome.'


def test_orffinder_valid_orfs(orf_list):
    global total_orfs, valid_orfs

    valid_orfs = 0
    total_orfs = 999999

    check_init()

    num_student = len(student_orfs)
    student_set = set(student_orfs)

    passed = (num_student > 0)

    assert passed, 'List is empty.'

    passed = (num_student == len(student_set))

    assert passed, 'List contains duplicates.'

    valid_orfs = len((orf_list & student_set))
    total_orfs = len((orf_list | student_set))

    passed = (valid_orfs > 0)

    assert passed, 'List does not contain any valid orfs.'


# The following tests will always fail if the ORF list contains duplicates!
def test_orffinder_valid_orfs_25():
    passed = (valid_orfs >= ceil(0.25 * total_orfs))
    assert passed, 'Less than 25% valid ORFs.'


def test_orffinder_valid_orfs_50():
    passed = (valid_orfs >= ceil(0.50 * total_orfs))
    assert passed, 'Less than 50% valid ORFs.'


def test_orffinder_valid_orfs_75():
    passed = (valid_orfs >= ceil(0.75 * total_orfs))
    assert passed, 'Less than 75% valid ORFs.'


def test_orffinder_valid_orfs_95():
    passed = (valid_orfs >= ceil(0.95 * total_orfs))
    assert passed, 'Less than 95% valid ORFs.'


def test_orffinder_valid_orfs_100_x1():
    passed = (valid_orfs == total_orfs)
    assert passed, 'Less than 100% valid ORFs.'


def test_orffinder_valid_orfs_100_x2():
    passed = (valid_orfs == total_orfs)
    assert passed, 'Less than 100% valid ORFs.'


def test_orffinder_valid_orfs_100_x3():
    passed = (valid_orfs == total_orfs)
    assert passed, 'Less than 100% valid ORFs.'
