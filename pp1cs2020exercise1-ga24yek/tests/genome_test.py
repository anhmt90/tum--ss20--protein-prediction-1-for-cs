import json
import pytest

from math import isclose
from genome import Genome
from pathlib import Path


@pytest.fixture(scope="module")
def json_data():
    test_json = 'genome_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def sequence(json_data):
    return json_data['sequence']


@pytest.fixture(scope="module")
def gc_content(json_data):
    return json_data['gc_content']


@pytest.fixture(scope="module")
def codon_dist(json_data):
    return json_data['codon_dist']


@pytest.fixture(scope="module")
def amino_acid_dist(json_data):
    return json_data['amino_acid_dist']


genome = None
passed_codon_test = False
passed_amino_acid_test = False


def check_init():
    assert genome is not None, 'Genome initialization failed.'


def test_genome_get_gc_content(sequence, gc_content):
    global genome

    passed = True

    try:
        genome = Genome(sequence)
    except Exception:
        genome = None
        passed = False

    assert passed, 'Error while creating Genome.'

    check_init()

    try:
        gc = genome.get_gc_content()
    except Exception:
        passed = False

    assert passed, 'Error in Genome.get_gc_content().'

    try:
        passed = isclose(gc, gc_content, rel_tol=0, abs_tol=1e-5)
    except Exception:
        passed = False

    assert passed, 'Incorrect GC content.'


def compare_dicts(obj_a, obj_b):
    if type(obj_a) != type(dict()) or type(obj_b) != type(dict()):
        try:
            return isclose(obj_a, obj_b, rel_tol=0, abs_tol=1e-5)
        except Exception:
            return False

    keys_a = set(obj_a.keys())
    keys_b = set(obj_b.keys())

    if keys_a != keys_b:
        return False

    return all([compare_dicts(obj_a[key], obj_b[key]) for key in keys_a])


def test_genome_get_codon_dist(codon_dist):
    global passed_codon_test

    passed_codon_test = False

    check_init()

    passed = True

    try:
        cd = genome.get_codon_dist()
    except Exception:
        passed = False
    finally:
        passed_codon_test = False

    assert passed, 'Error in Genome.get_codon_dist().'

    passed = compare_dicts(codon_dist, cd)

    assert passed, 'Incorrect codon distribution.'

    passed_codon_test = True


def test_genome_get_codon_dist_bonus():
    assert passed_codon_test, 'Did not pass test_genome_get_codon_dist().'


def test_genome_get_amino_acid_dist(amino_acid_dist):
    global passed_amino_acid_test

    passed_amino_acid_test = False

    check_init()

    passed = True

    try:
        aad = genome.get_amino_acid_dist()
    except Exception:
        passed = False
    finally:
        passed_amino_acid_test = False

    assert passed, 'Error in Genome.get_amino_acid_dist().'

    passed = compare_dicts(amino_acid_dist, aad)

    assert passed, 'Incorrect amino acid distribution.'

    passed_amino_acid_test = True


def test_genome_get_amino_acid_dist_bonus():
    assert passed_amino_acid_test, 'Did not pass test_genome_get_amino_acid_dist().'
