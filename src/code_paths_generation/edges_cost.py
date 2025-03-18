import numpy as np


def t_func_int(xij, qij, tij=1, alpha=0.15, beta=4):

    return tij*(xij + (alpha*xij*pow(xij*1.0/qij, beta)/(beta + 1)))

def t_func(xij, qij, tij=1, alpha=0.15, beta=4):

    return tij*(1.0 + (alpha*pow(xij*1.0/qij, beta)))

