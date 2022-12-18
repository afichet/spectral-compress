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

#include "Util.h"

#include <vector>
#include <string>
#include <cstring>

void Util::split_extension(const char* filename, std::string& base_str, std::string& extension_str)
{
    std::vector<char> base, extension;

    const size_t filename_len = std::strlen(filename);

    int len_base = filename_len;

    while (len_base >= 0 && filename[len_base] != '.') {
        --len_base;
    }

    if (len_base == -1) {
        // No extension
        base.resize(filename_len + 1);
        std::memcpy(base.data(), filename, sizeof(char) * filename_len);
        base.back() = 0;

        extension.resize(1);
        extension[0] = 0;
    } else {
        base.resize(len_base + 1);
        std::memcpy(base.data(), filename, sizeof(char) * len_base);
        base[len_base] = 0;

        extension.resize(filename_len - len_base + 1);
        std::memcpy(extension.data(), &filename[len_base], sizeof(char) * extension.size());
    }

    base_str = std::string(base.data());
    extension_str = std::string(extension.data());
}


