
/**
 *
 * author: Arthur Dufay @ inria.fr
 * Copyright INRIA 2017
 *
 * Modifications:
 *
 * Alban Fichet @ institutoptique.fr
 * April 2020
 * Copyright CNRS 2020
 *
 **/

#pragma once

#include <string>
#include <functional>
#include <vector>
#include <map>

class ENVISpectralImage
{
public:
    /**
     * Order for saving the data in the .raw file
     */
    enum InterleavedMode
    {
        IM_BIL = 0,
        IM_BIP,
        IM_BSQ,
        IM_UNDEF
    };

    /**
     * Type of data stored in the .raw file
     */
    enum DataType
    {
        DT_UINT_8 = 1,
        DT_INT_16 = 2,
        DT_INT_32 = 3,
        DT_FLOAT = 4,
        DT_DOUBLE = 5,
        DT_CPLX_FLOAT = 6,
        DT_CPLX_DOUBLE = 9,
        DT_UINT_16 = 12,
        DT_UINT_32 = 13,
        DT_INT_64 = 14,
        DT_UINT_64 = 15,
        DT_UNDEF
    };

    /**
     * @class Header
     * Store information about the specific attributes of an ENVI file
     * contained in the header part (.hdr)
     */
    class Header
    {
    public:
        Header(
            size_t _width = 0,
            size_t _height = 0,
            size_t _bands = 0,
            InterleavedMode _interleaved_mode = IM_BSQ,
            std::string const &_comments = "",
            DataType _data_type = DT_FLOAT)
            : width(_width)
            , height(_height)
            , bands(_bands)
            , interleave_mode(_interleaved_mode)
            , comments(_comments)
            , data_type(_data_type)
        {
        }

        Header(Header const &cpy)
            : width(cpy.width)
            , height(cpy.height)
            , bands(cpy.bands)
            , interleave_mode(cpy.interleave_mode)
            , comments(cpy.comments)
            , data_type(cpy.data_type)
        {
        }

        size_t width;
        size_t height;
        size_t bands;
        InterleavedMode interleave_mode;
        std::string comments;
        DataType data_type;
    };

    /**
     * @brief Creates a new ENVISpectralImage from a stored file
     *
     * We assume the header file and the image data share the same name
     * exept for the extension: <name>.hdr and <name>.raw
     *
     * @param file_path Either a .raw or a .hdr file, the second file will
     *                  be used using the same filename with the other
     *                  extension.
     * @param progress  A callback to notify the loading progress in
     *                  percents.
     * @throws IMAGE_LOAD_SAVE_FLAGS if unsucessfull
     */
    ENVISpectralImage(
        const std::string &file_path,
        std::function<void(int)> progress = [](int) {});

    virtual ~ENVISpectralImage();

    size_t width() const { return _header.width; }
    size_t height() const { return _header.height; }

    /**
     * @brief Get the pair of filename for a given ENVI file
     *
     * ENVI format is using two files: one header (.hdr) and one file
     * containing the image data (.raw).
     * This method gives the two corresponding filenames given one filename
     *
     * @param file_path  One of the filename (either .hdr or .raw).
     * @param names      Result:
     *                      - name.first is the header (.hdr).
     *                      - name.second is the data (.raw).
     *
     * @returns true if the provided file_path is either ending with .hdr
     *          or .raw. false otherwise.
     */
    static bool getFilenames(std::string const &file_path, std::pair<std::string, std::string> &names);

    std::vector<float>& getImage() { return _image_data; }
    const std::vector<float>& wavelengths() const { return _wavelengths; }

private:
    //---------------------------------------------------------------------
    // Read routines
    //---------------------------------------------------------------------

    /**
     * @brief Loads in ENVI file format
     *
     * In this method, we assume the header file and the image data have the
     * same name exept for the extension: <name>.hdr and <name>.raw
     *
     * @param file_path Either a .raw or a .hdr file, the second file will
     *                  be used using the same filename with the other
     *                  extension.
     * @param progress  A callback to notify the loading progress in
     *                  percents.
     *
     * @return The load state. true if successfull.
     */
    bool load(
        const std::string &file_path,
        std::function<void(int)> progress = [](int) {});

    /**
     * @brief Loads in ENVI file format
     *
     * @param file_path_header The header (.hdr) file.
     * @param file_path_data   The data (.raw) file.
     * @param progress  A callback to notify the loading progress in
     *                  percents.
     *
     * @return The load state. true if successfull.
     */
    bool load(
        const std::string &file_path_header,
        const std::string &file_path_image_data,
        std::function<void(int)> progress = [](int) {});

    /**
     * @brief Reads the header from an input stream.
     *
     * @param is Input stream to the ENVI file header.
     * @returns The load state. true if successfull.
     */
    bool readHeader(std::istream &is);

    /**
     * @brief Reads data from an output stream.
     *
     * @param is        Input stream to the ENVI file data.
     * @param progress  A callback to notify the loading progress in
     *                  percents.
     *
     * @returns The load state. true if successfull.
     */
    template <typename T>
    bool readImageData(
        std::istream &is,
        std::function<void(int)> progress = [](int) {});

    // FIELDS
    Header _header;
    std::vector<float> _wavelengths;
    std::vector<float> _image_data;
};
