"""Unit tests of SONAT"""

ORDERED_MODULES = ["fcore", "pack", "stack", "misc", "plot", "render",
    "ens", "obs", "arm", "cui"]


def print_test_order(mode="module"):

    assert mode in ('module', 'file', 'prefix')

    if mode=='prefix':
        print ' '.join([('test_'+m) for m in ORDERED_MODULES])

    elif mode=='module':
        print ' '.join([('sonat.test.' + p) for p in ORDERED_MODULES])
    else:
        print ' '.join([(p + '.py') for p in ORDERED_MODULES])


