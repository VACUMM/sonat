"""Test script for module :mod:`sonat.my`"""

import os
from util import THIS_DIR


from sonat.my import load_user_code_file
from sonat.obs import get_obs_platform

def test_my_load_user_code_file():

    load_user_code_file(os.path.join(THIS_DIR, 'mysonat.py'))



if __name__=='__main__':
    test_my_load_user_code_file()




