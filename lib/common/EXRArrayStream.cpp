/**
 * Copyright 2022 Alban Fichet
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

#include "EXRArrayStream.h"

#include <cstring>

EXRArrayStream::EXRArrayStream()
: Imf::OStream("mem")
, Imf::IStream("mem")
, _pos(0)
{}


EXRArrayStream::EXRArrayStream(const std::vector<uint8_t> data)
: Imf::OStream("mem")
, Imf::IStream("mem")
, _data(data)
, _pos(0)
{}


void EXRArrayStream::write(const char c[/*n*/], int n) {
    const uint64_t remaining_bytes = _data.size() - _pos;

    if ((int)remaining_bytes < n) {
        _data.resize(_data.size() + n - remaining_bytes);
    }

    std::memcpy(&_data[_pos], c, n);

    _pos += n;
}


bool EXRArrayStream::read(char c[/*n*/], int n) {
    const uint64_t remaining_bytes = _data.size() - _pos;

    if ((int)remaining_bytes < n) {
        throw std::exception();
    }

    std::memcpy(c, &_data[_pos], n);

    _pos += n;

    return _pos == _data.size();
}


uint64_t EXRArrayStream::tellp() {
    return _pos;
}


uint64_t EXRArrayStream::tellg() {
    return _pos;
}


void EXRArrayStream::seekp(uint64_t pos) {
    _pos = pos;
}


void EXRArrayStream::seekg(uint64_t pos) {
    _pos = pos;
}


const std::vector<uint8_t>& EXRArrayStream::data() const {
    return _data;
}


size_t EXRArrayStream::size() const {
    return _data.size();
}
