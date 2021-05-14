import exe0_subseq

import pytest

def test_subseq(data_case_ex0_subseq ):
    result_inc, result_dec = exe0_subseq.get_longest_subsequences(data_case_ex0_subseq["input"])

    test_result = result_inc == data_case_ex0_subseq["output"][0] and result_dec == data_case_ex0_subseq["output"][1]
    assert test_result