import ctypes
import os

from . import consts

# def get_dll(path, *args):
#     try:
#         return ctypes.CDLL(path, *args)
#     except OSError:
#         print(
#             textwrap.dedent(
#                 f"""
#         Failed to import mcl shared library from:

#         {DIR_FOR_LINKER}

#         Please set your mcl installation dir path to MCL_PATH and run again

#         export MCL_PATH=<path to mcl library>

#         """
#             )
#         )
#         sys.exit(1)


# with change_cwd(DIR_FOR_LINKER):
#     system = platform.system()
#     if system == "Darwin":
#         mclbls12_384 = get_dll("lib/libmclbn384_256.dylib")
#     elif system == "Linux":
#         get_dll("lib/libmcl.so", ctypes.RTLD_GLOBAL)
#         mclbls12_384 = get_dll("lib/libmclbn384_256.so")
#     else:
#         raise RuntimeError(f"Unsupported OS {system}")

import platform
if platform.system()=="Windows":
    lib_mcl_path=os.path.join(os.path.dirname(__file__),"mcl_384_256.dll")
    mcl=ctypes.WinDLL(lib_mcl_path, ctypes.RTLD_GLOBAL)
    mclbls12_384=mcl
elif platform.system()=="Linux":
    lib_mcl_path=os.path.join(os.path.dirname(__file__),"lishe384_256.so")
    mcl=ctypes.CDLL(lib_mcl_path, ctypes.RTLD_GLOBAL)
    mclbls12_384=mcl
else:
    raise RuntimeError(platform.system())
ret = mcl.mclBn_init(consts.BLS12_381, consts.MCLBN_COMPILED_TIME_VAR)

if ret:
    raise RuntimeError(f"mclbls12_384 ret {ret}")
