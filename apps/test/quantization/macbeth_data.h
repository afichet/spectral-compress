/**
 * Copyright (c) 2020 - 2021
 * Alban Fichet, Romain Pacanowski, Alexander Wilkie
 * Institut d'Optique Graduate School, CNRS - Universite de Bordeaux,
 * Inria, Charles University
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  * Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *  * Neither the name of Institut d'Optique Graduate School, CNRS -
 * Universite de Bordeaux, Inria, Charles University nor the names of
 * its contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#pragma once

// Data from http://www.babelcolor.com/colorchecker-2.htm

const float macbeth_wavelengths[36]
  = {380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490,
     500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610,
     620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730};
const float macbeth_patches[24][36] = {
  {0.055, 0.058, 0.061, 0.062, 0.062, 0.062, 0.062, 0.062, 0.062,
   0.062, 0.062, 0.063, 0.065, 0.070, 0.076, 0.079, 0.081, 0.084,
   0.091, 0.103, 0.119, 0.134, 0.143, 0.147, 0.151, 0.158, 0.168,
   0.179, 0.188, 0.190, 0.186, 0.181, 0.182, 0.187, 0.196, 0.209},
  {0.117, 0.143, 0.175, 0.191, 0.196, 0.199, 0.204, 0.213, 0.228,
   0.251, 0.280, 0.309, 0.329, 0.333, 0.315, 0.286, 0.273, 0.276,
   0.277, 0.289, 0.339, 0.420, 0.488, 0.525, 0.546, 0.562, 0.578,
   0.595, 0.612, 0.625, 0.638, 0.656, 0.678, 0.700, 0.717, 0.734},
  {0.130, 0.177, 0.251, 0.306, 0.324, 0.330, 0.333, 0.331, 0.323,
   0.311, 0.298, 0.285, 0.269, 0.250, 0.231, 0.214, 0.199, 0.185,
   0.169, 0.157, 0.149, 0.145, 0.142, 0.141, 0.141, 0.141, 0.143,
   0.147, 0.152, 0.154, 0.150, 0.144, 0.136, 0.132, 0.135, 0.147},
  {0.051, 0.054, 0.056, 0.057, 0.058, 0.059, 0.060, 0.061, 0.062,
   0.063, 0.065, 0.067, 0.075, 0.101, 0.145, 0.178, 0.184, 0.170,
   0.149, 0.133, 0.122, 0.115, 0.109, 0.105, 0.104, 0.106, 0.109,
   0.112, 0.114, 0.114, 0.112, 0.112, 0.115, 0.120, 0.125, 0.130},
  {0.144, 0.198, 0.294, 0.375, 0.408, 0.421, 0.426, 0.426, 0.419,
   0.403, 0.379, 0.346, 0.311, 0.281, 0.254, 0.229, 0.214, 0.208,
   0.202, 0.194, 0.193, 0.200, 0.214, 0.230, 0.241, 0.254, 0.279,
   0.313, 0.348, 0.366, 0.366, 0.359, 0.358, 0.365, 0.377, 0.398},
  {0.136, 0.179, 0.247, 0.297, 0.320, 0.337, 0.355, 0.381, 0.419,
   0.466, 0.510, 0.546, 0.567, 0.574, 0.569, 0.551, 0.524, 0.488,
   0.445, 0.400, 0.350, 0.299, 0.252, 0.221, 0.204, 0.196, 0.191,
   0.188, 0.191, 0.199, 0.212, 0.223, 0.232, 0.233, 0.229, 0.229},
  {0.054, 0.054, 0.053, 0.054, 0.054, 0.055, 0.055, 0.055, 0.056,
   0.057, 0.058, 0.061, 0.068, 0.089, 0.125, 0.154, 0.174, 0.199,
   0.248, 0.335, 0.444, 0.538, 0.587, 0.595, 0.591, 0.587, 0.584,
   0.584, 0.590, 0.603, 0.620, 0.639, 0.655, 0.663, 0.663, 0.667},
  {0.122, 0.164, 0.229, 0.286, 0.327, 0.361, 0.388, 0.400, 0.392,
   0.362, 0.316, 0.260, 0.209, 0.168, 0.138, 0.117, 0.104, 0.096,
   0.090, 0.086, 0.084, 0.084, 0.084, 0.084, 0.084, 0.085, 0.090,
   0.098, 0.109, 0.123, 0.143, 0.169, 0.205, 0.244, 0.287, 0.332},
  {0.096, 0.115, 0.131, 0.135, 0.133, 0.132, 0.130, 0.128, 0.125,
   0.120, 0.115, 0.110, 0.105, 0.100, 0.095, 0.093, 0.092, 0.093,
   0.096, 0.108, 0.156, 0.265, 0.399, 0.500, 0.556, 0.579, 0.588,
   0.591, 0.593, 0.594, 0.598, 0.602, 0.607, 0.609, 0.609, 0.610},
  {0.092, 0.116, 0.146, 0.169, 0.178, 0.173, 0.158, 0.139, 0.119,
   0.101, 0.087, 0.075, 0.066, 0.060, 0.056, 0.053, 0.051, 0.051,
   0.052, 0.052, 0.051, 0.052, 0.058, 0.073, 0.096, 0.119, 0.141,
   0.166, 0.194, 0.227, 0.265, 0.309, 0.355, 0.396, 0.436, 0.478},
  {0.061, 0.061, 0.062, 0.063, 0.064, 0.066, 0.069, 0.075, 0.085,
   0.105, 0.139, 0.192, 0.271, 0.376, 0.476, 0.531, 0.549, 0.546,
   0.528, 0.504, 0.471, 0.428, 0.381, 0.347, 0.327, 0.318, 0.312,
   0.310, 0.314, 0.327, 0.345, 0.363, 0.376, 0.381, 0.378, 0.379},
  {0.063, 0.063, 0.063, 0.064, 0.064, 0.064, 0.065, 0.066, 0.067,
   0.068, 0.071, 0.076, 0.087, 0.125, 0.206, 0.305, 0.383, 0.431,
   0.469, 0.518, 0.568, 0.607, 0.628, 0.637, 0.640, 0.642, 0.645,
   0.648, 0.651, 0.653, 0.657, 0.664, 0.673, 0.680, 0.684, 0.688},
  {0.066, 0.079, 0.102, 0.146, 0.200, 0.244, 0.282, 0.309, 0.308,
   0.278, 0.231, 0.178, 0.130, 0.094, 0.070, 0.054, 0.046, 0.042,
   0.039, 0.038, 0.038, 0.038, 0.038, 0.039, 0.039, 0.040, 0.041,
   0.042, 0.044, 0.045, 0.046, 0.046, 0.048, 0.052, 0.057, 0.065},
  {0.052, 0.053, 0.054, 0.055, 0.057, 0.059, 0.061, 0.066, 0.075,
   0.093, 0.125, 0.178, 0.246, 0.307, 0.337, 0.334, 0.317, 0.293,
   0.262, 0.230, 0.198, 0.165, 0.135, 0.115, 0.104, 0.098, 0.094,
   0.092, 0.093, 0.097, 0.102, 0.108, 0.113, 0.115, 0.114, 0.114},
  {0.050, 0.049, 0.048, 0.047, 0.047, 0.047, 0.047, 0.047, 0.046,
   0.045, 0.044, 0.044, 0.045, 0.046, 0.047, 0.048, 0.049, 0.050,
   0.054, 0.060, 0.072, 0.104, 0.178, 0.312, 0.467, 0.581, 0.644,
   0.675, 0.690, 0.698, 0.706, 0.715, 0.724, 0.730, 0.734, 0.738},
  {0.058, 0.054, 0.052, 0.052, 0.053, 0.054, 0.056, 0.059, 0.067,
   0.081, 0.107, 0.152, 0.225, 0.336, 0.462, 0.559, 0.616, 0.650,
   0.672, 0.694, 0.710, 0.723, 0.731, 0.739, 0.746, 0.752, 0.758,
   0.764, 0.769, 0.771, 0.776, 0.782, 0.790, 0.796, 0.799, 0.804},
  {0.145, 0.195, 0.283, 0.346, 0.362, 0.354, 0.334, 0.306, 0.276,
   0.248, 0.218, 0.190, 0.168, 0.149, 0.127, 0.107, 0.100, 0.102,
   0.104, 0.109, 0.137, 0.200, 0.290, 0.400, 0.516, 0.615, 0.687,
   0.732, 0.760, 0.774, 0.783, 0.793, 0.803, 0.812, 0.817, 0.825},
  {0.108, 0.141, 0.192, 0.236, 0.261, 0.286, 0.317, 0.353, 0.390,
   0.426, 0.446, 0.444, 0.423, 0.385, 0.337, 0.283, 0.231, 0.185,
   0.146, 0.118, 0.101, 0.090, 0.082, 0.076, 0.074, 0.073, 0.073,
   0.074, 0.076, 0.077, 0.076, 0.075, 0.073, 0.072, 0.074, 0.079},
  {0.189, 0.255, 0.423, 0.660, 0.811, 0.862, 0.877, 0.884, 0.891,
   0.896, 0.899, 0.904, 0.907, 0.909, 0.911, 0.910, 0.911, 0.914,
   0.913, 0.916, 0.915, 0.916, 0.914, 0.915, 0.918, 0.919, 0.921,
   0.923, 0.924, 0.922, 0.922, 0.925, 0.927, 0.930, 0.930, 0.933},
  {0.171, 0.232, 0.365, 0.507, 0.567, 0.583, 0.588, 0.590, 0.591,
   0.590, 0.588, 0.588, 0.589, 0.589, 0.591, 0.590, 0.590, 0.590,
   0.589, 0.591, 0.590, 0.590, 0.587, 0.585, 0.583, 0.580, 0.578,
   0.576, 0.574, 0.572, 0.571, 0.569, 0.568, 0.568, 0.566, 0.566},
  {0.144, 0.192, 0.272, 0.331, 0.350, 0.357, 0.361, 0.363, 0.363,
   0.361, 0.359, 0.358, 0.358, 0.359, 0.360, 0.360, 0.361, 0.361,
   0.360, 0.362, 0.362, 0.361, 0.359, 0.358, 0.355, 0.352, 0.350,
   0.348, 0.345, 0.343, 0.340, 0.338, 0.335, 0.334, 0.332, 0.331},
  {0.105, 0.131, 0.163, 0.180, 0.186, 0.190, 0.193, 0.194, 0.194,
   0.192, 0.191, 0.191, 0.191, 0.192, 0.192, 0.192, 0.192, 0.192,
   0.192, 0.193, 0.192, 0.192, 0.191, 0.189, 0.188, 0.186, 0.184,
   0.182, 0.181, 0.179, 0.178, 0.176, 0.174, 0.173, 0.172, 0.171},
  {0.068, 0.077, 0.084, 0.087, 0.089, 0.090, 0.092, 0.092, 0.091,
   0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090, 0.090,
   0.090, 0.090, 0.090, 0.089, 0.089, 0.088, 0.087, 0.086, 0.086,
   0.085, 0.084, 0.084, 0.083, 0.083, 0.082, 0.081, 0.081, 0.081},
  {0.031, 0.032, 0.032, 0.033, 0.033, 0.033, 0.033, 0.033, 0.032,
   0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032,
   0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032,
   0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.032, 0.033}};