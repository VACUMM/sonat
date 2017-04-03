"""Unit tests of SONAT"""


ORDER = "test_fcore test_pack test_stack test_misc test_plot test_render test_ens test_obs test_arm test_cui"

def print_test_order(mode="module"):

    assert mode in ('module', 'file', 'prefix')

    if mode=='prefix':
        print ORDER

    elif mode=='module':
        print ' '.join([('sonat.test.' + p) for p in ORDER.split(' ')])
    else:
        print ' '.join([(p + '.py') for p in ORDER.split(' ')])


