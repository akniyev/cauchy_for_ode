//
//  laguerre.c
//  laguerre
//
//  Created by Hasan Akniyev on 27.11.2019.
//  Copyright © 2019 Hasan Akniyev. All rights reserved.
//

#include "laguerre.h"

double dabs(double x) {
    if (x < 0) return -x;
    return x;
}

double laguerre_iterative(int k, double alpha, double x) {
    if (k <= 0) {
        return 1;
    }
    
    double result = -x + alpha + 1;
    
    double minus_2 = 1;
    double minus_1 = result;
    
    for (int i = 2; i <= k; i++) {
        double a = (2.0 * i - 1 + alpha - x) / i;
        double b = (i + alpha - 1.0) / i;
        
        double current = a * minus_1 - b * minus_2;
        
        minus_2 = minus_1;
        minus_1 = current;
    }
    
    result = minus_1;
    
    return result;
}

double root_of_laguerre(int k, double alpha, double a, double b, double epsilon) {
    double mid_point = (a + b) / 2.0;
    double fa = laguerre_iterative(k, alpha, a);
    double fc = laguerre_iterative(k, alpha, mid_point);
    
    if (dabs(b-a) < epsilon || dabs(fc) < epsilon) {
        return mid_point;
    }
    
    if (fa * fc < 0) {
        return root_of_laguerre(k, alpha, a, mid_point, epsilon);
    } else {
        return root_of_laguerre(k, alpha, mid_point, b, epsilon);
    }
}
