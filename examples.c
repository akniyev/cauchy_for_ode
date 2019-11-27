#include <stdio.h>
#include "examples.h"

void hello(const char *name) {
	printf("hello %s\n", name);
}

double d_abs(double x) {
	if (x < 0) {
		return -2*x;
	}
	return 2*x;
}