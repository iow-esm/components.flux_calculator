import sys
import os
try:
    module_dir = str(sys.argv[1]) 
except:
    module_dir = os.getcwd()

import cffi
ffibuilder = cffi.FFI()


header = """
extern int set(char *, int32_t *);
extern int get(char *, int32_t *);
extern int call_function(char *, char *);
"""

module = f"""
from pyfort import ffi
import numpy as np
import importlib
import os
import sys

# Create the dictionary mapping ctypes to np dtypes.
ctype2dtype = {{}}

# Integer types
for prefix in ('int', 'uint'):
    for log_bytes in range(4):
        ctype = '%s%d_t' % (prefix, 8 * (2**log_bytes))
        dtype = '%s%d' % (prefix[0], 2**log_bytes)
        #print( ctype )
        #print( dtype )
        ctype2dtype[ctype] = np.dtype(dtype)

# Floating point types
ctype2dtype['float'] = np.dtype('f4')
ctype2dtype['double'] = np.dtype('f8')


def asarray(ffi, ptr, shape, **kwargs):
    length = np.prod(shape)
    # Get the canonical C type of the elements of ptr as a string.
    T = ffi.getctype(ffi.typeof(ptr).item)
    #print( T )
    #print( ffi.sizeof( T ) )

    if T not in ctype2dtype:
        raise RuntimeError("Cannot create an array for element type: %s" % T)

    a = np.frombuffer(ffi.buffer(ptr, length * ffi.sizeof(T)), ctype2dtype[T]).reshape(shape, **kwargs)
    return a

STATE = {{}}

@ffi.def_extern(error=1)
def get(key, c_ptr, shape=(1,)):
    \"\"\"Copy the numpy array stored in STATE[key] to a pointer\"\"\"

    key = ffi.string(key).decode("UTF-8")

    # wrap pointer in numpy array
    fortran_arr = asarray(ffi, c_ptr, shape)
    
    # update the numpy array in place
    fortran_arr[:] = STATE[key]
  
    return 0
    
@ffi.def_extern(error=1)
def set(key, c_ptr, shape=(1,)):
    \"\"\"Call python\"\"\"
    key = ffi.string(key).decode("UTF-8")
    STATE[key] = asarray(ffi, c_ptr, shape).copy()

    return 0

	
@ffi.def_extern(error=1)
def call_function(module_name, function_name):
    #pwd = os.getcwd()
    #print(pwd)
    #sys.path.append(pwd)
    sys.path.append("{module_dir}")
    module_name = ffi.string(module_name).decode("UTF-8")
    function_name = ffi.string(function_name).decode("UTF-8")

    mod = importlib.import_module(module_name)
    
    # the function we want to call
    fun = getattr(mod, function_name)
    
    # call the function
    # this function can edit STATE inplace
    fun(STATE)

    return 0
"""

with open("pyfort.h", "w") as f:
    f.write(header)

ffibuilder.embedding_api(header)
ffibuilder.set_source("pyfort", r'''
    #include "pyfort.h"
''',)

ffibuilder.embedding_init_code(module)
ffibuilder.emit_c_code("pyfort.c")
ffibuilder.compile(target="libpyfort.so", verbose=True)
