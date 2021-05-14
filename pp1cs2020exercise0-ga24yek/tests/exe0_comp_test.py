import exe0_comp

import pytest


def test_complementary(data_case_ex0_comp):
    result = exe0_comp.complementary(data_case_ex0_comp["input"]).upper()
    test_result = result == data_case_ex0_comp["output"]
    assert test_result