import numpy as np
import pandas as pd
import dimod
import neal


def get_feasible(A, b, sampler, samples=20):

    AA = np.dot(A.T, A)
    h = -2.0*np.dot(b.T, A)
    Q = AA + np.diag(h[0])
    offset = np.dot(b.T, b) + 0.0

    # Define Binary Quadratic Model
    bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)

    response = sampler.sample(bqm, num_reads=samples)

    response = response.aggregate()

    filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]

    feas_sols = response.record.sample[filter_idx]

    return feas_sols


def get_feasible_encode(A, b, E, L, sampler, samples=20):

    Q_I = np.dot(A.T, A)
    h = 2.0*np.dot(np.dot(L.T, Q_I) - np.dot(b.T, A), E)
    Q = np.dot(E.T, np.dot(Q_I, E)) + np.diag(h[0])
    offset = np.dot(L.T, np.dot(Q_I, L)) + np.dot(b.T, b) - 2.0*np.dot(b.T, np.dot(A, L))

    # Define Binary Quadratic Model
    bqm = dimod.BinaryQuadraticModel.from_numpy_matrix(mat=Q, offset=offset)

    response = sampler.sample(bqm, num_reads=samples)

    response = response.aggregate()

    filter_idx = [i for i, e in enumerate(response.record.energy) if e == 0.0]

    feas_sols_X = response.record.sample[filter_idx]

    print(L.shape, E.shape, feas_sols_X.shape)
    feas_sols = (L + np.dot(E, feas_sols_X.T)).T

    return feas_sols_X, feas_sols, Q_I, h, Q, offset


