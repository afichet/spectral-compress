/**
 * Copyright 2022 - 2023 Alban Fichet
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer in the documentation and/or other materials provided
 * with the distribution.
 *
 *   3. Neither the name of the copyright holder nor the names of its
 * contributors may be used to endorse or promote products derived
 * from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 * STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "cmd_constraints.h"


std::string FrameDistanceConstraint::description() const
{
    return "Sets the distance level for lossy compression";
}


std::string FrameDistanceConstraint::shortID() const
{
    return "0..15";
}


bool FrameDistanceConstraint::check(const float &value) const
{
    return (value >= 0.f) && (value <= 15.f);
}


// ----------------------------------------------------------------------------


std::string CompressionEffortConstraint::description() const
{
    return "Sets the compression effort for JXL.";
}


std::string CompressionEffortConstraint::shortID() const
{
    return "1..9";
}


bool CompressionEffortConstraint::check(const int &value) const
{
    return (value >= 1) && (value <= 9);
}


// ----------------------------------------------------------------------------


std::string DownsamplingFactorConstraint::description() const
{
    return "Sets the spatial downsampling ratio to apply to AC components.";
}


std::string DownsamplingFactorConstraint::shortID() const
{
    return "1, 2, 4, or 8";
}


bool DownsamplingFactorConstraint::check(const int &value) const
{
    return value == 1
        || value == 2
        || value == 4
        || value == 8;
}
