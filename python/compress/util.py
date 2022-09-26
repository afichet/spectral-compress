# Copyright (c) 2019, Christoph Peters
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Karlsruhe Institute of Technology nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import numpy as np

def levinsons_algorithm(first_column: np.array) -> np.array:
    """Implements Levinson's algorithm to solve a special system of linear 
       equations. The matrix has the given first column and is a Hermitian 
       Toeplitz matrix, i.e. it has constant diagonals. The right-hand side 
       vector is the canonical basis vector (1,0,...,0).
      \return A vector of same shape as first_column.
      \note This algorithm provides an efficient way to compute the vector q in 
            the paper, which is called evaluation_polynomial in this program."""
    solution = np.zeros_like(first_column)
    solution[0] = 1.0 / first_column[0]

    for j in range(1, first_column.shape[0]):
        dot_product = np.dot(solution[0:j], first_column[j:0:-1])
        solution[0:j + 1] = (solution[0:j + 1] - dot_product * np.conj(solution[j::-1])) / (1 - np.abs(dot_product) ** 2)

    return solution


def levinsons_algorithm_bias(first_column: np.array) -> np.array:
    """Implements Levinson's algorithm to solve a special system of linear 
       equations. The matrix has the given first column and is a Hermitian 
       Toeplitz matrix, i.e. it has constant diagonals. The right-hand side 
       vector is the canonical basis vector (1,0,...,0).
      \param CorrectionBias Pass a small positive constant (e.g. 1.0e-4) to 
             enable biasing. The vector FirstColumn is considered invalid if 
             the resulting Toeplitz matrix has a negative eigenvalue. In this 
             case, entries of FirstColumn will be modified on the fly to make 
             the matrix positive definite.
      \return A vector of same shape as FirstColumn.
      \note This algorithm provides an efficient way to compute the vector q in 
            the paper, which is called EvaluationPolynomial in this program."""
    solution = np.zeros_like(first_column)
    solution[0] = 1.0 / first_column[0]
    epsilon = 1e-4

    for j in range(1, first_column.shape[0]):
        partial_dot_product = np.dot(solution[1:j], first_column[j-1:0:-1])
        dot_product = solution[0] * first_column[j] + partial_dot_product

        if np.abs(dot_product) >= 1.0:
            dot_product *= (1.0 - epsilon) / np.abs(dot_product)
            first_column[j] = (dot_product - partial_dot_product) / solution[0]
            epsilon=1.0

        solution[0:j + 1] = (solution[0:j+1] - dot_product * np.conj(solution[j::-1])) / (1 - np.abs(dot_product)**2)

    return solution


def get_levinson_dots(first_column: np.array) -> np.array:
    """
    :return: A vector holding the dot products that arise in each iteration of
        Levinson's algorithm.
    """
    solution = np.zeros_like(first_column)
    solution[0] = 1.0 / first_column[0]

    dots = np.zeros_like(first_column)

    for j in range(1, first_column.shape[0]):
        dots[j] = np.dot(solution[0:j], first_column[j:0:-1])
        solution[0:j + 1] = (solution[0:j + 1] - dots[j] * np.conj(solution[j::-1])) / (1 - np.abs(dots[j]) ** 2)
    
    return dots


def run_levinson_from_dots(diagonal_entry: complex, dots: np.array) -> np.array:
    """
    Runs Levinson's algorithm using the given intermediate dot products and the
    diagonal entry. Outputs the corresponding input vector and the result of
    Levinson's algorithm.
    """
    solution = np.zeros_like(dots)
    first_column = np.zeros_like(dots)

    solution[0]     = 1.0 / diagonal_entry
    first_column[0] = diagonal_entry

    for j in range(1, dots.shape[0]):
        radius = 1.0 / solution[0].real
        center = -radius * np.dot(solution[1:j], first_column[j - 1:0:-1])
        first_column[j] = center + radius * dots[j]
        solution[0:j + 1] = (solution[0:j + 1] - dots[j] * np.conj(solution[j::-1])) / (1 - np.abs(dots[j]) ** 2)
    
    return first_column



def normalize(vec: np.array):
    mins = np.min(vec, axis=0)
    maxs = np.max(vec, axis=0)

    vec_rescaled = np.zeros_like(vec)

    vec_rescaled[:, 0] = vec[:, 0]

    for i in range(1, vec_rescaled.shape[1]):
        vec_rescaled[:, i] = (vec[:, i] - mins[i]) / (maxs[i] - mins[i])

    return vec_rescaled, mins, maxs


def denormalize(vec_rescaled: np.array, mins: np.array, maxs: np.array) -> np.array:
    vec = np.zeros_like(vec_rescaled)

    vec[:, 0] = vec_rescaled[:, 0]

    for i in range(1, vec.shape[1]):
        vec[:, i] = vec_rescaled[:, i] * (maxs[i] - mins[i]) + mins[i]

    return vec
