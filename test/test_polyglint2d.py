from setiptah.polyglint2d import *
import itertools
import numpy as np

import pytest


def test_polyglint2d():

    # triangle
    Ab = np.array([
        [1, 1, 1],
        [-1, 0, 0],
        [0, -1, 0],
    ])

    A = Ab[:, :-1]
    b = Ab[:, -1]

    # unit density
    c = np.array([0, 0])
    d = 1

    assert integrate(c, d, A, b) == .5


def test_vertices_of_triangle():

    Ab = np.array([
        [1, 1, 1],
        [-1, 0, 0],
        [0, -1, 0],
    ])

    A = Ab[:, :-1]
    b = Ab[:, -1]

    vertices = enumerate_vertices_2d(A, b)

    expected = [
        (0, 0),
        (1, 0),
        (0, 1),
    ]

    def score_match(observed):
        return sum(
            np.linalg.norm(o - np.array(e), ord=2)
            for o, e in zip(observed, expected)
        )

    best = min(
        score_match(observed)
        for observed in itertools.permutations(vertices)
    )

    assert best < .00001
