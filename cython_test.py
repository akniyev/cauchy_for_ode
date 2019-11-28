from laguerre_functions import *


k = 100
alpha = 0

start_time = time.time()
roots = find_all_roots_of_laguerre(k, alpha)
end_time = time.time()

print(roots)
print(end_time - start_time)
