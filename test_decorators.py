import time


def cached(func):
    memo = {}

    def wrapper(*args):
        if args in memo:
            return memo[args]
        result = func(*args)
        memo[args] = result
        return result

    return wrapper


@cached
def f(a):
    time.sleep(1)
    return a ** 2


# @cached
# def g(a):
#     time.sleep(1)
#     return a ** 3


print(f(1))
print(f(2))
print(f(3))

print(f(1))
print(f(2))
print(f(3))

# print(g(1))
# print(g(2))
# print(g(3))
