
import sys
sys.path.insert(0, '..')


from sonat.obs import NcObsPlatform


class MyObsPlatform(NcObsPlatform):
    platform_type = 'myplatform'
    pass
