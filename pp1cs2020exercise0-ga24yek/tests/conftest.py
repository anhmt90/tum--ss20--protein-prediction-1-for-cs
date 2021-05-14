import pytest
import os
import sys
import json


def json_data(file):
    with open(os.path.join(os.path.dirname(__file__), file)) as json_file:
        json_data = json.load(json_file)
    return json_data

def pytest_generate_tests(metafunc):
    if "data_case_ex0_comp" in metafunc.fixturenames:
        case_data = json_data('exe0_comp_test.json')
        cases = [case_data[i] for i in case_data]
        metafunc.parametrize("data_case_ex0_comp", cases)

    if "data_case_ex0_sub" in metafunc.fixturenames:
        case_data = json_data('exe0_sub_test.json')
        cases = [case_data[i] for i in case_data]
        metafunc.parametrize("data_case_ex0_sub", cases)
    
    if "data_case_ex0_subseq" in metafunc.fixturenames:
        case_data = json_data('exe0_subseq_test.json')
        cases = [case_data[i] for i in case_data]
        metafunc.parametrize("data_case_ex0_subseq", cases)