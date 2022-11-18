#include "ENVISpectralImage.h"

#include <fstream>
#include <stdexcept>
#include <regex>
#include <cassert>
#include <limits>
#include <cmath>
#include <sstream>

ENVISpectralImage::ENVISpectralImage(
    const std::string &file_path,
    std::function<void(int)> progress)
    : _header()
{
    if (!load(file_path, progress)) {
        throw std::runtime_error("Could not load ENVI image");
    }
}


ENVISpectralImage::~ENVISpectralImage() {}


//---------------------------------------------------------------------
// Read routines
//---------------------------------------------------------------------

bool ENVISpectralImage::load(const std::string &file_path, std::function<void(int)> progress)
{
    std::pair<std::string, std::string> filenames;

    if (getFilenames(file_path, filenames)) {
        return load(filenames.first, filenames.second, progress);
    } else {
        return false;
    }
}


bool ENVISpectralImage::load(
    std::string const &file_path_header,
    std::string const &file_path_image_data,
    std::function<void(int)> progress)
{
    bool status;

    std::ifstream header_ifs(file_path_header, std::ifstream::in);

    if (!header_ifs.is_open()) {
        return false;
    }

    status = readHeader(header_ifs);
    header_ifs.close();

    // If an error occured while writing the header, we stop there
    if (!status) {
        return false;
    }

    std::ifstream data_ifs(file_path_image_data, std::ios::in | std::ios::binary | std::ios::ate);

    if (!data_ifs.is_open()) {
        return false;
    }

    switch (_header.data_type) {
    case DT_UINT_8:  status = readImageData<uint8_t> (data_ifs, progress); break;
    case DT_UINT_16: status = readImageData<uint16_t>(data_ifs, progress); break;
    case DT_UINT_32: status = readImageData<uint32_t>(data_ifs, progress); break;
    case DT_UINT_64: status = readImageData<uint64_t>(data_ifs, progress); break;
    case DT_INT_16:  status = readImageData<int16_t> (data_ifs, progress); break;
    case DT_INT_32:  status = readImageData<int32_t> (data_ifs, progress); break;
    case DT_INT_64:  status = readImageData<int64_t> (data_ifs, progress); break;
    case DT_FLOAT:   status = readImageData<float>   (data_ifs, progress); break;
    case DT_DOUBLE:  status = readImageData<double>  (data_ifs, progress); break;
    case DT_CPLX_FLOAT:
    case DT_CPLX_DOUBLE:
    case DT_UNDEF:
        data_ifs.close();
        return false;
    }

    data_ifs.close();

    return status;
}


bool ENVISpectralImage::getFilenames(
    const std::string &file_path,
    std::pair<std::string, std::string> &names)
{
    std::string extension = file_path.substr(file_path.size() - 3, file_path.size());
    std::string base_name = file_path.substr(0, file_path.size() - 4);

    if (extension == "hdr" || extension == "HDR") {
        // User specified the header, we guess the data
        names.first = file_path;
        names.second = base_name + ".raw";

        return true;
    } else if (extension == "raw" || extension == "RAW") {
        // User specified the data, we guess the header
        names.first = base_name + ".hdr";
        names.second = file_path;

        return true;
    }

    return false;
}


bool ENVISpectralImage::readHeader(std::istream &is)
{
    std::string line;

    std::function<bool(std::string const &, std::string const &)> check_param = [](
                                                                                    std::string const &line,
                                                                                    std::string const &param_name)
    {
        return line.substr(0, param_name.size()) == param_name;
    };

    std::function<std::string(std::string const &)> get_param_value = [](std::string const &line)
    {
        std::size_t param_position = line.find_first_of('=');
        std::string param_value;

        if (param_position != std::string::npos)
        {
            param_value = line.substr(param_position + 1);
            param_value = std::regex_replace(param_value, std::regex("^ +"), "$1");
        }
        return param_value;
    };

    while (std::getline(is, line)) {
        // remove leading spaces
        line = std::regex_replace(line, std::regex("^ +"), "$1");
        std::string param_value = "";

        // read description
        if (check_param(line, "description")) {
            if (line.find_first_of('{') != std::string::npos)
            {
                line = line.substr(line.find_first_of('{') + 1);
                while (line.find_first_of('}') == std::string::npos && !is.eof()) {
                    _header.comments += line + "\n";
                    std::getline(is, line);
                }
                std::replace(line.begin(), line.end(), '}', ' ');
                _header.comments += line;
            }
        }
        // file type
        else if (check_param(line, "file type")) {
            param_value = get_param_value(line);

            if (!check_param(param_value, "ENVI")) {
                return false;
            }
        }
        // interleave method
        else if (check_param(line, "interleave")) {
            param_value = get_param_value(line);

            if (check_param(param_value, "bil")) {
                _header.interleave_mode = IM_BIL;
            }
            else if (check_param(param_value, "bsq")) {
                _header.interleave_mode = IM_BSQ;
            }
            else if (check_param(param_value, "bip")) {
                _header.interleave_mode = IM_BIP;
            }
            else {
                _header.interleave_mode = IM_UNDEF;
                return false;
            }
        }
        //
        else if (check_param(line, "samples")) {
            param_value = get_param_value(line);
            _header.width = stoi(param_value);
        }
        //
        else if (check_param(line, "bands")) {
            param_value = get_param_value(line);
            _header.bands = stoi(param_value);
        }
        //
        else if (check_param(line, "lines")) {
            param_value = get_param_value(line);
            _header.height = stoi(param_value);
        }
        else if (check_param(line, "data type")) {
            param_value = get_param_value(line);
            int data_type = stoi(param_value);

            switch (data_type) {
            case DT_UINT_8:
            case DT_UINT_16:
            case DT_UINT_32:
            case DT_UINT_64:
            case DT_INT_16:
            case DT_INT_32:
            case DT_INT_64:
            case DT_FLOAT:
            case DT_DOUBLE:
            case DT_CPLX_FLOAT:
            case DT_CPLX_DOUBLE:
                _header.data_type = static_cast<DataType>(data_type);
                break;
            default:
                _header.data_type = DT_UNDEF;
                return false;
            }
        }
        else if (check_param(line, "byte order")) {
            param_value = get_param_value(line);
            if (stoi(param_value) != 0) {
                return false;
            }
        }
        // read wavelength
        else if (check_param(line, "Wavelength") || check_param(line, "wavelength")) {
            if (line.find_first_of('{') != std::string::npos) {
                bool end_wave = false;

                while (!end_wave && !is.eof()) {
                    std::getline(is, line);

                    if (line.find_first_of('}') != std::string::npos) {
                        end_wave = true;
                        line.erase(std::remove(line.begin(), line.end(), '}'), line.end());
                    }

                    // std::vector<std::string> tokens;
                    // mrf::util::StringParsing::tokenize(line, tokens, ",");

                    auto checkSpace = [](unsigned char const c) {
                        return std::isspace(c);
                    };

                    std::string str;
                    std::stringstream stream(line);

                    while (std::getline(stream, str, ',')) {
                        str.erase(std::remove_if(str.begin(), str.end(), checkSpace), str.end());
                        
                        if (str.size() > 0) {
                            float wavelength_f = std::round(stof(str));
                            _wavelengths.push_back(static_cast<int>(wavelength_f));
                        }
                    }
                }
            }
        }
    }

    return true;
}


template <typename T>
bool ENVISpectralImage::readImageData(std::istream &is, std::function<void(int)> progress)
{
    _image_data.resize(width() * height() * _wavelengths.size());

    const size_t buffer_size = _image_data.size();
    const size_t n_bands = _wavelengths.size();

    std::unique_ptr<T> buffer(new T[buffer_size]);
    size_t read_offset = 0;

    is.seekg(0, is.end);
    size_t stream_size = is.tellg();
    is.seekg(0, is.beg);

    if (stream_size < buffer_size * sizeof(T)) {
        return false;
    }

    is.seekg(0, std::ios::beg);
    is.read((char *)(buffer.get()), buffer_size * sizeof(T));

    if (!is) {
        return false;
    }

    switch (_header.interleave_mode) {
    case IM_BIL:
        for (uint i = 0; i < height(); i++) {
            for (uint slice = 0; slice < n_bands; slice++) {
                for (uint j = 0; j < width(); j++) {
                    assert(read_offset < buffer_size);

                    if (std::is_integral<T>::value) {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]) / static_cast<float>(std::numeric_limits<T>::max());
                    }
                    else {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]);
                    }
                }
            }

            progress(int(100.f * float(i) / float(height() - 1)));
        }
        break;

    case IM_BIP:
        for (uint i = 0; i < height(); i++) {
            for (uint j = 0; j < width(); j++) {
                for (uint slice = 0; slice < n_bands; slice++) {
                    assert(read_offset < buffer_size);

                    if (std::is_integral<T>::value)
                    {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]) / static_cast<float>(std::numeric_limits<T>::max());
                    }
                    else
                    {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]);
                    }
                }
            }

            progress(int(100.f * float(i) / float(height() - 1)));
        }
        break;

    case IM_BSQ:
        for (uint slice = 0; slice < n_bands; slice++) {
            for (uint i = 0; i < height(); i++) {
                for (uint j = 0; j < width(); j++) {
                    assert(read_offset < buffer_size);

                    if (std::is_integral<T>::value) {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]) / static_cast<float>(std::numeric_limits<T>::max());
                    }
                    else {
                        _image_data[n_bands * (i * width() + j) + slice] =
                            static_cast<float>(buffer.get()[read_offset++]);
                    }
                }
            }

            progress(int(100.f * float(slice) / float(n_bands - 1)));
        }
        break;

    default:
        return false;
    }

    return true;
}
