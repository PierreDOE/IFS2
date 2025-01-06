# -*- coding: utf-8 -*-
"""
tools to do something general
"""
import numpy as np
from viv.colored_messages import color_s, set_msg

def test_complex_numbers(name, u):
    """
    :param name: a name of the table
    :param u: a table of dimension (n, 2)
    :return: the first indexes where the table is complex
    """
    set_msg("Research : complex number  for %s" % name)
    n, m = u.shape
    for j in range(m):
        for i in range(n):
            if np.iscomplex(u[i, 1]):
                print("\t component : %i, complex number at index :" % j, color_s(str(i), "red"))
                print("\t \t values: ", u[i, :])
                break
    print("\t" + color_s("max", "red"), np.amax(np.abs(u), axis=0), np.argmax(np.abs(u), axis=0))
    print("\t" + color_s("min", "blue"), np.amin(np.abs(u), axis=0), np.argmin(abs(np.abs(u)), axis=0))

def test_real_parts(name, u):
    """
    :param name: a name of the table
    :param u: a table of dimension (n, 2)
    :return: the first indexes where the real part of the table is zero
    """
    set_msg("Reseach: a real part == 0 for %s " % name)
    n, m = u.shape
    for j in range(m):
        for i in range(n):
            if abs(np.real(u[i, j])) <= 1e-10:
                print("\t component : %i, pure imaginary number at index :" % j, color_s(str(i), "red"))
                print("\t \t values: ", u[i, :])
                break

    print("\t" + color_s("max", "red"), np.amax(np.real(u), axis=0), np.argmax(np.real(u), axis=0))
    print("\t" + color_s("min", "blue"), np.amin(np.real(u), axis=0), np.argmin(abs(np.real(u)), axis=0))

def test_imaginary_parts(name, u):
    """
    test if complex and delete imaginary part
    a complex array  N x 2
    :param name: a name of the table
    :param u: a table of dimension (n, 2)
    :return: the first indexes where the imaginary part of the table is zero
    """
    set_msg("Research: imaginary part == 0 for %s " %name)
    n, m = u.shape
    for j in range(m):
        for i in range(n):
            if abs(np.imag(u[i, j])) <= 1e-10:
                print("\t component : %i, pure imaginary number at index :" % j, color_s(str(i), "red"))
                print("\t \t values: ", u[i, :])
                break

    print("\t" + color_s("max", "red"), np.amax(np.imag(u), axis=0), np.argmax(np.imag(u), axis=0))
    print("\t" + color_s("min", "blue"), np.amin(np.imag(u), axis=0), np.argmin(abs(np.imag(u)), axis=0))


def test_values(name, u):
    """
    test values of an array u
    :param name:
    :param u:
    :return:
    """
    test_complex_numbers(name, u)
    test_real_parts(name, u)
    test_imaginary_parts(name, u)

def test_array_is_nan(name, u):
    """
    find nan value in a array
    a complex array  N x 2
    :param name: a name of the table
    :param u: a table of dimension (n, 2)
    :return: the first indexes where the imaginary part of the table is zero
    """
    set_msg("Research:  nan values for %s " % name)
    n, m = u.shape
    s = np.isnan(u)
    for j in range(m):
        for i in range(n):
            if s[i, j]:
                print("\t component : %i, is nan at index :" % j, color_s(str(i), "red"))
                print("\t \t values: ", u[i, :])


