import json
import pytest
import numpy as np

from pathlib import Path
from .matrices import MATRICES
from local_alignment import LocalAlignment


@pytest.fixture(scope="module")
def json_data():
    test_json = 'local_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def null(json_data):
    return json_data['null']


@pytest.fixture(scope="module")
def short(json_data):
    return json_data['short']


null_la = None
short_la = None


def create_LA(string_1, string_2, gap_penalty, matrix):
    passed = True

    try:
        la = LocalAlignment(string_1, string_2, gap_penalty, matrix)
    except Exception:
        passed = False

    assert passed, 'Error while creating LocalAlignment.'

    assert la is not None, 'LocalAlignment initialization failed.'

    return la


def test_has_alignment(null, short):
    global null_la, short_la

    null_la = create_LA(*null['strings'],
                        null['gap_penalty'],
                        MATRICES[null['matrix']])

    short_la = create_LA(*short['strings'],
                         short['gap_penalty'],
                         MATRICES[short['matrix']])

    passed = True

    try:
        has_alignment = null_la.has_alignment()
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.has_alignment().'

    passed = (not has_alignment)

    assert passed, 'Incorrect value for has_alignment(); null strings.'

    try:
        has_alignment = short_la.has_alignment()
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.has_alignment().'

    passed = (has_alignment)

    assert passed, 'Incorrect value for has_alignment(); small strings.'


def test_get_alignment_on_small_strings(short):
    assert short_la is not None, 'LocalAlignment initialization failed.'

    passed = True

    try:
        alignment = short_la.get_alignment()
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.get_alignment().'

    passed = (tuple(short['alignment']) == alignment)

    assert passed, 'Incorrect alignment (small strings).'


def test_get_alignment_on_null_strings(null):
    assert null_la is not None, 'LocalAlignment initialization failed.'

    passed = True

    try:
        alignment = null_la.get_alignment()
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.get_alignment().'

    passed = (tuple(null['alignment']) == alignment)

    assert passed, 'Incorrect alignment (null strings).'


def test_is_residue_aligned_on_first_string(short):
    assert short_la is not None, 'LocalAlignment initialization failed.'

    [s_1, p_1, res_1], [s_2, p_2, res_2] = short["residue_aligned_on_first"]

    passed = True

    try:
        is_aligned_1 = short_la.is_residue_aligned(s_1, p_1)
        is_aligned_2 = short_la.is_residue_aligned(s_2, p_2)
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.is_residue_aligned().'

    passed = ((res_1 == is_aligned_1) and (res_2 == is_aligned_2))

    assert passed, 'Incorrect value for is_residue_aligned(); first string.'


def test_is_residue_aligned_on_second_string(short):
    assert short_la is not None, 'LocalAlignment initialization failed.'

    [s_1, p_1, res_1], [s_2, p_2, res_2] = short["residue_aligned_on_second"]

    passed = True

    try:
        is_aligned_1 = short_la.is_residue_aligned(s_1, p_1)
        is_aligned_2 = short_la.is_residue_aligned(s_2, p_2)
    except Exception:
        passed = False

    assert passed, 'Error in LocalAlignment.is_residue_aligned().'

    passed = ((res_1 == is_aligned_1) and (res_2 == is_aligned_2))

    assert passed, 'Incorrect value for is_residue_aligned(); second string.'
