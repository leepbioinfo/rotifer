# For now, we use this as part of our script's template
import os
import sys
_add_path = [ os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib"),
        os.path.join(os.path.dirname(os.path.dirname(__file__)), "lib", "python" + str(sys.version_info.major) + "." + str(sys.version_info.minor), "site-packages")
        ]
for _d in _add_path:
    if os.path.exists(_d):
        sys.path.insert(0,_d)
import rotifer.core.cli as corecli
