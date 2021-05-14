import exe0_sub
import pytest


def test_subsearch(data_case_ex0_sub, ):
    result = exe0_sub.get_sub_position(*data_case_ex0_sub["input"])
    test_result = result == data_case_ex0_sub["output"]
    assert test_result