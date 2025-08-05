import ctypes
from ..hook import mcl as mcl_lib
from .. import builder
from . import base
from .G1 import G1
from .G2 import G2
from .Fp import Fp

mcl_lib.mclBnGT_isEqual.argtypes=[
    ctypes.POINTER(base.Structure),
    ctypes.POINTER(base.Structure)
]
mcl_lib.mclBnGT_isEqual.restype=int

mcl_lib.mclBnGT_mul.argtypes=[
    ctypes.POINTER(base.Structure),
    ctypes.POINTER(base.Structure),
    ctypes.POINTER(base.Structure),
]
mcl_lib.mclBnGT_mul.restype=None

mcl_lib.mclBnGT_add.argtypes=[
    ctypes.POINTER(base.Structure),
    ctypes.POINTER(base.Structure),
    ctypes.POINTER(base.Structure),
]

mcl_lib.mclBnGT_add.restype=None

@builder.provide_methods(
    builder.method("__invert__").using(builder.buildTwoOp).with_args("inv"),
    builder.method("pairing").using(builder.buildPairing).with_args(G1, G2),
    builder.method("getStr"),
)


class GT(base.Structure):
    _fields_ = [
        ("d", (Fp * 12)),
    ]

    def __eq__(self, gt2: "GT") -> bool:
        return mcl_lib.mclBnGT_isEqual(self,gt2)==1

    def __mul__(self,gt2:"GT") -> "GT":
        ret=GT()
        # print("before",ret)
        mcl_lib.mclBnGT_mul(ret,self,gt2)
        # print("after",ret)
        return ret

    # def __add__(self,gt2:"GT") -> "GT":
    #     ret=GT()
    #     # print("before",ret)
    #     mcl_lib.mclBnGT_add(ret,self,gt2)
    #     # print("after",ret)
    #     return ret
