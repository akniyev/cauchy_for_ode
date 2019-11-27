import ctypes

lib = ctypes.cdll.LoadLibrary('./laguerre_cpp_lib/liblaguerre.dylib')

print(lib)

print(lib.__Z3sayNSt3__112basic_stringIcNS_11char_traitsIcEENS_9allocatorIcEEEE())
