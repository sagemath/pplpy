# content of test_sample.py
import os
def func(x):
    return x + 1

def test_answer():
    print(os.getcwd())
    print(os.listdir())
    assert func(3) == 6

def test_one():
    assert 1==1