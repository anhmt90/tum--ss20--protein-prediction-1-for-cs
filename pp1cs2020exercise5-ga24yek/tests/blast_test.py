import json
import pytest
import numpy as np

from time import time as ttime
from blast import BlastDB, Blast
from pathlib import Path


@pytest.fixture(scope="module")
def json_data():
    test_json = 'blast_test.json'
    relative_path = Path(__file__).parent

    with Path(relative_path, test_json).open('r') as json_file:
        json_data = json.load(json_file)

    return json_data


@pytest.fixture(scope="module")
def db_sequences(json_data):
    return json_data['db_sequences']


@pytest.fixture(scope="module")
def db_stats(json_data):
    return tuple(json_data['db_stats'])


@pytest.fixture(scope="module")
def db_seqs_for_word(json_data):
    return json_data['db_seqs_for_word']


@pytest.fixture(scope="module")
def sub_matrix(json_data):
    return np.array(json_data['sub_matrix'], dtype=np.int64)


@pytest.fixture(scope="module")
def query_seq(json_data):
    return json_data['query_seq']


@pytest.fixture(scope="module")
def query_word(json_data):
    return json_data['query_word']


@pytest.fixture(scope="module")
def query_pssm(json_data):
    return np.array(json_data['query_pssm'], dtype=np.int64)


@pytest.fixture(scope="module")
def blast_words(json_data):
    return json_data['blast_words']


@pytest.fixture(scope="module")
def blast_words_pssm(json_data):
    return json_data['blast_words_pssm']


def table_list_tuple(data):
    for key, value in data.items():
        data[key] = [tuple(x) for x in value]
    return data


@pytest.fixture(scope="module")
def blast_hsp_one_hit(json_data):
    return table_list_tuple(json_data['blast_hsp_one_hit'])


@pytest.fixture(scope="module")
def blast_hsp_two_hit(json_data):
    return table_list_tuple(json_data['blast_hsp_two_hit'])


@pytest.fixture(scope="module")
def blast_hsp_one_hit_pssm(json_data):
    return table_list_tuple(json_data['blast_hsp_one_hit_pssm'])


@pytest.fixture(scope="module")
def blast_hsp_two_hit_pssm(json_data):
    return table_list_tuple(json_data['blast_hsp_two_hit_pssm'])


blast = None
blast_db = None


def check_init():
    assert blast is not None, 'Blast initialization failed.'
    assert blast_db is not None, 'BlastDB initialization failed.'


@pytest.mark.timeout(30)
def test_blast_db_get_db_stats(db_sequences, db_stats):
    global blast_db

    passed = True

    try:
        blast_db = BlastDB()
    except Exception:
        passed = False
        blast_db = None

    assert passed, 'Error while creating BlastDB.'

    assert blast_db is not None, 'BlastDB initialization failed.'

    try:
        for s in db_sequences:
            blast_db.add_sequence(s)
    except Exception:
        passed = False

    assert passed, 'Error in BlastDB.add_sequence().'

    try:
        stats = blast_db.get_db_stats()
    except Exception:
        passed = False

    assert passed, 'Error in BlastDB.get_db_stats().'

    passed = (db_stats == stats)

    assert passed, 'Incorrect BlastDB statistics.'


@pytest.mark.timeout(5)
def test_blast_db_get_sequences(query_word, db_seqs_for_word):
    assert blast_db is not None, 'BlastDB initialization failed.'

    passed = True

    try:
        seqs_for_word = blast_db.get_sequences(query_word)
    except Exception:
        passed = False

    assert passed, 'Error in BlastDB.get_sequences().'

    try:
        passed_1 = (len(db_seqs_for_word) == len(seqs_for_word))
        passed_2 = (set(db_seqs_for_word) == set(seqs_for_word))
    except Exception:
        passed = False

    assert passed, 'Error while comparing BlastDB.get_sequences() output.'

    passed = (passed_1 and passed_2)

    assert passed, 'Incorrect sequences returned.'


@pytest.mark.timeout(5)
def test_blast_get_words(sub_matrix, query_seq, blast_words):
    global blast

    passed = True

    try:
        blast = Blast(sub_matrix)
    except Exception:
        blast = None
        passed = False

    assert passed, 'Error while creating Blast.'

    assert blast is not None, 'Blast initialization failed.'

    try:
        words = blast.get_words(sequence=query_seq, T=13)
    except Exception:
        passed = False

    assert passed, 'Error in Blast.get_words(sequence).'

    try:
        passed_1 = (len(blast_words) == len(words))
        passed_2 = (set(blast_words) == set(words))
    except Exception:
        passed = False

    assert passed, 'Error while comparing Blast.get_words(sequence) output.'

    passed = (passed_1 and passed_2)

    assert passed, 'Incorrect words returned for sequence.'


@pytest.mark.timeout(5)
def test_blast_get_words_with_pssm(query_pssm, blast_words_pssm):
    assert blast is not None, 'Blast initialization failed.'

    passed = True

    try:
        words = blast.get_words(pssm=query_pssm, T=11)
    except Exception:
        passed = False

    assert passed, 'Error in Blast.get_words(pssm).'

    try:
        passed_1 = (len(blast_words_pssm) == len(words))
        passed_2 = (set(blast_words_pssm) == set(words))
    except Exception:
        passed = False

    assert passed, 'Error while comparing Blast.get_words(pssm) output.'

    passed = (passed_1 and passed_2)

    assert passed, 'Incorrect words returned for PSSM.'


def compare_blast_results(blast_results, results, hint):
    passed = True

    try:
        passed_1 = (len(blast_results) == len(results))
        passed_2 = (set(blast_results) == set(results))
    except Exception:
        passed = False

    assert passed, f'Error while comparing Blast results ({hint}).'

    passed = (passed_1 and passed_2)

    assert passed, f'Incorrect target sequences returned ({hint}).'

    for target, hsp_list in results.items():
        blast_hsp_list = blast_results[target]

        try:
            passed_1 = (len(blast_hsp_list) == len(hsp_list))
            passed_2 = (set(blast_hsp_list) == set(hsp_list))
        except Exception:
            passed = False

        assert passed, f'Error while comparing Blast HSP list ({hint}).'

        passed = (passed_1 and passed_2)

        assert passed, f'Incorrect HSPs returned ({hint}).'


run_time_one_hit = 9999
run_time_two_hit = 9999
run_time_one_hit_pssm = 9999
run_time_two_hit_pssm = 9999


@pytest.mark.timeout(60)
def test_blast_search_one_hit(query_seq, blast_hsp_one_hit):
    global run_time_one_hit

    run_time_one_hit = 9999

    check_init()

    passed = True

    run_time_t1 = ttime()

    try:
        results = blast.search_one_hit(blast_db,
                                       query=query_seq,
                                       T=13,
                                       X=5,
                                       S=30)
    except Exception:
        passed = False
    finally:
        run_time_one_hit = 9999

    run_time_t2 = ttime()

    assert passed, 'Error in Blast.search_one_hit(sequence)'

    compare_blast_results(blast_hsp_one_hit, results, 'one-hit, sequence')

    run_time_one_hit = int(run_time_t2 - run_time_t1)


def test_blast_search_one_hit_30s():
    assert run_time_one_hit >= 0, 'Invalid runtime.'
    assert run_time_one_hit <= 30, f'Too slow (one-hit, sequence): {run_time_one_hit}s'


def test_blast_search_one_hit_15s():
    assert run_time_one_hit >= 0, 'Invalid runtime.'
    assert run_time_one_hit <= 15, f'Too slow (one-hit, sequence): {run_time_one_hit}s'


@pytest.mark.timeout(60)
def test_blast_search_two_hit(query_seq, blast_hsp_two_hit):
    global run_time_two_hit

    run_time_two_hit = 9999

    check_init()

    passed = True

    run_time_t1 = ttime()

    try:
        results = blast.search_two_hit(blast_db,
                                       query=query_seq,
                                       T=11,
                                       X=5,
                                       S=30,
                                       A=40)
    except Exception:
        passed = False
    finally:
        run_time_two_hit = 9999

    run_time_t2 = ttime()

    assert passed, 'Error in Blast.search_two_hit(sequence)'

    compare_blast_results(blast_hsp_two_hit, results, 'two-hit, sequence')

    run_time_two_hit = int(run_time_t2 - run_time_t1)


def test_blast_search_two_hit_20s():
    assert run_time_two_hit >= 0, 'Invalid runtime.'
    assert run_time_two_hit <= 20, f'Too slow (two-hit, sequence): {run_time_two_hit}s'


def test_blast_search_two_hit_10s():
    assert run_time_two_hit >= 0, 'Invalid runtime.'
    assert run_time_two_hit <= 10, f'Too slow (two-hit, sequence): {run_time_two_hit}s'


@pytest.mark.timeout(60)
def test_blast_search_one_hit_with_pssm(query_pssm, blast_hsp_one_hit_pssm):
    global run_time_one_hit_pssm

    run_time_one_hit_pssm = 9999

    check_init()

    passed = True

    run_time_t1 = ttime()

    try:
        results = blast.search_one_hit(blast_db,
                                       pssm=query_pssm,
                                       T=13,
                                       X=5,
                                       S=30)
    except Exception:
        passed = False
    finally:
        run_time_one_hit_pssm = 9999

    run_time_t2 = ttime()

    assert passed, 'Error in Blast.search_one_hit(pssm)'

    compare_blast_results(blast_hsp_one_hit_pssm, results, 'one-hit, pssm')

    run_time_one_hit_pssm = int(run_time_t2 - run_time_t1)


def test_blast_search_one_hit_with_pssm_30s():
    assert run_time_one_hit_pssm >= 0, 'Invalid runtime.'
    assert run_time_one_hit_pssm <= 30, f'Too slow (one-hit, pssm): {run_time_one_hit_pssm}s'


@pytest.mark.timeout(60)
def test_blast_search_two_hit_with_pssm(query_pssm, blast_hsp_two_hit_pssm):
    global run_time_two_hit_pssm

    run_time_two_hit_pssm = 9999

    check_init()

    passed = True

    run_time_t1 = ttime()

    try:
        results = blast.search_two_hit(blast_db,
                                       pssm=query_pssm,
                                       T=11,
                                       X=5,
                                       S=30,
                                       A=40)
    except Exception:
        passed = False
    finally:
        run_time_two_hit_pssm = 9999

    run_time_t2 = ttime()

    assert passed, 'Error in Blast.search_two_hit(pssm)'

    compare_blast_results(blast_hsp_two_hit_pssm, results, 'two-hit, pssm')

    run_time_two_hit_pssm = int(run_time_t2 - run_time_t1)


def test_blast_search_two_hit_with_pssm_20s():
    assert run_time_two_hit_pssm >= 0, 'Invalid runtime.'
    assert run_time_two_hit_pssm <= 20, f'Too slow (two-hit, pssm): {run_time_two_hit_pssm}s'


def test_blast_search_all_15s():
    run_time_max = max(run_time_one_hit,
                       run_time_two_hit,
                       run_time_one_hit_pssm,
                       run_time_two_hit_pssm)

    assert run_time_max >= 0, 'Invalid runtime.'
    assert run_time_max <= 15, f'Too slow (max runtime): {run_time_max}s'
