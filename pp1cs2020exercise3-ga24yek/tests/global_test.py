import json
import pytest
import numpy as np

from pathlib import Path
from .matrices import MATRICES
from global_alignment import GlobalAlignment


@pytest.fixture(scope="module")
def json_data():
    test_json = 'global_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def small(json_data):
    return json_data['small']


@pytest.fixture(scope="module")
def large(json_data):
    return json_data['large']


@pytest.fixture(scope="module")
def matching(json_data):
    return json_data['matching']


@pytest.fixture(scope="module")
def mismatching(json_data):
    return json_data['mismatching']


@pytest.fixture(scope="module")
def indel(json_data):
    return json_data['indel']


@pytest.fixture(scope="module")
def all_changes(json_data):
    return json_data['all_changes']


small_ga = None
large_ga = None


def create_GA(string_1, string_2, gap_penalty, matrix):
    passed = True

    try:
        ga = GlobalAlignment(string_1, string_2, gap_penalty, matrix)
    except Exception:
        passed = False

    assert passed, 'Error while creating GlobalAlignment.'

    assert ga is not None, 'GlobalAlignment initialization failed.'

    return ga


def test_get_best_score_on_small_strings(small):
    global small_ga

    small_ga = create_GA(*small['strings'],
                         small['gap_penalty'],
                         MATRICES[small['matrix']])

    passed = True

    try:
        score = small_ga.get_best_score()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_best_score().'

    passed = (small['best_score'] == score)

    assert passed, 'Incorrect best score (small strings).'


def test_get_best_score_on_large_strings(large):
    global large_ga

    large_ga = create_GA(*large['strings'],
                         large['gap_penalty'],
                         MATRICES[large['matrix']])

    passed = True

    try:
        score = large_ga.get_best_score()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_best_score().'

    passed = (large['best_score'] == score)

    assert passed, 'Incorrect best score (large strings).'


def test_get_number_of_alignments_on_small_strings(small):
    assert small_ga is not None, 'GlobalAlignment initialization failed.'

    passed = True

    try:
        num_alignments = small_ga.get_number_of_alignments()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_number_of_alignments().'

    passed = (small['number_of_alignments'] == num_alignments)

    assert passed, 'Incorrect number of alignments (small strings).'


def test_get_number_of_alignments_on_large_strings(large):
    assert large_ga is not None, 'GlobalAlignment initialization failed.'

    passed = True

    try:
        num_alignments = large_ga.get_number_of_alignments()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_number_of_alignments().'

    passed = (large['number_of_alignments'] == num_alignments)

    assert passed, 'Incorrect number of alignments (large strings).'


def test_get_alignments_on_small_strings(small):
    assert small_ga is not None, 'GlobalAlignment initialization failed.'

    passed = True

    try:
        alignments = small_ga.get_alignments()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_alignments().'

    passed = ({tuple(x) for x in small['alignments']} == set(alignments))

    assert passed, 'Incorrect alignments (small strings).'


def test_get_alignments_on_large_strings(large):
    assert large_ga is not None, 'GlobalAlignment initialization failed.'

    passed = True

    try:
        alignments = large_ga.get_alignments()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_alignments().'

    passed = ({tuple(x) for x in large['alignments']} == set(alignments))

    assert passed, 'Incorrect alignments (large strings).'


def test_score_matrix_on_matching_strings(matching):
    ga = create_GA(*matching['strings'],
                   matching['gap_penalty'],
                   MATRICES[matching['matrix']])

    passed = True

    try:
        score_matrix = ga.get_score_matrix()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_score_matrix().'

    try:
        passed = np.array_equal(np.array(score_matrix),
                                np.array(matching['score_matrix']))
    except Exception:
        passed = False

    assert passed, 'Incorrect score matrix (matching strings).'


def test_score_matrix_on_mismatching_strings(mismatching):
    ga = create_GA(*mismatching['strings'],
                   mismatching['gap_penalty'],
                   MATRICES[mismatching['matrix']])

    passed = True

    try:
        score_matrix = ga.get_score_matrix()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_score_matrix().'

    try:
        passed = np.array_equal(np.array(score_matrix),
                                np.array(mismatching['score_matrix']))
    except Exception:
        passed = False

    assert passed, 'Incorrect score matrix (mismatching strings).'


def test_score_matrix_on_indel_strings(indel):
    ga = create_GA(*indel['strings'],
                   indel['gap_penalty'],
                   MATRICES[indel['matrix']])

    passed = True

    try:
        score_matrix = ga.get_score_matrix()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_score_matrix().'

    try:
        passed = np.array_equal(np.array(score_matrix),
                                np.array(indel['score_matrix']))
    except Exception:
        passed = False

    assert passed, 'Incorrect score matrix (indel strings).'


def test_score_matrix_on_all_changes_strings(all_changes):
    ga = create_GA(*all_changes['strings'],
                   all_changes['gap_penalty'],
                   MATRICES[all_changes['matrix']])

    passed = True

    try:
        score_matrix = ga.get_score_matrix()
    except Exception:
        passed = False

    assert passed, 'Error in GlobalAlignment.get_score_matrix().'

    try:
        passed = np.array_equal(np.array(score_matrix),
                                np.array(all_changes['score_matrix']))
    except Exception:
        passed = False

    assert passed, 'Incorrect score matrix (all_changes strings).'
