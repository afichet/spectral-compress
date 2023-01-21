#include <iostream>
#include <vector>
#include <sstream>

#include <OpenEXR/ImfInputFile.h>
#include <OpenEXR/ImfOutputFile.h>
#include <OpenEXR/ImfChannelList.h>
#include <OpenEXR/ImfStringAttribute.h>
#include <OpenEXR/ImfFrameBuffer.h>
#include <OpenEXR/ImfHeader.h>

#include <OpenEXR/ImfAttribute.h>
#include <OpenEXR/ImfIO.h>


class ArrayStream : public Imf::OStream, public Imf::IStream {
public:
    ArrayStream()
    : Imf::OStream("mem")
    , Imf::IStream("mem")
    , _pos(0)
    {}


    ArrayStream(const std::vector<uint8_t> data)
    : Imf::OStream("mem")
    , Imf::IStream("mem")
    , _data(data)
    , _pos(0)
    {}


    virtual void write(const char c[/*n*/], int n) {
        const uint64_t remaining_bytes = _data.size() - _pos;

        if (remaining_bytes < n) {
            _data.resize(_data.size() + n - remaining_bytes);
        }

        std::memcpy(&_data[_pos], c, n);

        _pos += n;
    }


    virtual bool read(char c[/*n*/], int n) {
        const uint64_t remaining_bytes = _data.size() - _pos;

        if (remaining_bytes < n) {
            throw std::exception();
        }

        std::memcpy(c, &_data[_pos], n);

        _pos += n;

        return _pos == _data.size();
    }


    virtual uint64_t tellp () {
        return _pos;
    }


    virtual uint64_t tellg() {
        return _pos;
    }


    virtual void seekp(uint64_t pos) {
        _pos = pos;
    }


    virtual void seekg(uint64_t pos) {
        _pos = pos;
    }


    uint8_t* data() {
        return _data.data();
    }


    size_t size() const {
        return _data.size();
    }

private:
    std::vector<uint8_t> _data;
    uint64_t _pos;
};


int main(int argc, char* argv[]) {
    Imf::InputFile      exr_in(argv[1]);
    const Imf::Header&  exr_header     = exr_in.header();
    const Imath::Box2i& exr_datawindow = exr_header.dataWindow();

    uint32_t width  = exr_datawindow.max.x - exr_datawindow.min.x + 1;
    uint32_t height = exr_datawindow.max.y - exr_datawindow.min.y + 1;

    ArrayStream attr_stream;

    for (Imf::Header::ConstIterator it = exr_header.begin(); it != exr_header.end(); it++) {
        const char* attribute_name = it.name();
        const char* attribute_type = it.attribute().typeName();

        attr_stream.write(attribute_name, std::strlen(attribute_name) + 1);
        attr_stream.write(attribute_type, std::strlen(attribute_type) + 1);

        it.attribute().writeValueTo(attr_stream, 1);
    }

    // Saving metadata
    std::FILE *f = std::fopen("attr.dat", "wb");

    if (f) {
        std::fwrite(attr_stream.data(), 1, attr_stream.size(), f);
        std::fclose(f);
    }

    // Writing EXR
    attr_stream.seekg(0);

    Imf::Header       exr_header_out(width, height);
    Imf::ChannelList &exr_channels_out = exr_header_out.channels();
    Imf::FrameBuffer  exr_framebuffer_out;

    do {
        char c;
        std::stringstream attribute_name_stream;
        std::stringstream attribute_type_stream;

        do {
            attr_stream.read(&c, 1);
            attribute_name_stream << c;
        } while(c != 0);

        do {
            attr_stream.read(&c, 1);
            attribute_type_stream << c;
        } while (c != 0);

        const std::string attribute_name = attribute_name_stream.str();
        const std::string attribute_type = attribute_type_stream.str();

        Imf::Attribute* attr = Imf::Attribute::newAttribute(attribute_type.c_str());
        attr->readValueFrom(attr_stream, attr_stream.size() - attr_stream.tellg(), 1);

        if (std::strcmp(attribute_name.c_str(), "channels") != 0) {
            exr_header_out.insert(attribute_name, *attr);
        }

        delete attr;
    } while (!(attr_stream.tellg() == attr_stream.size()));

    std::vector<float> foo_fb(width * height);

    exr_channels_out.insert("Y", Imf::Channel(Imf::FLOAT));
    exr_framebuffer_out.insert("Y", Imf::Slice(Imf::FLOAT, (char*)foo_fb.data(), sizeof(float), width * sizeof(float)));

    Imf::OutputFile exr_out("Test.exr", exr_header_out);
    exr_out.setFrameBuffer(exr_framebuffer_out);
    exr_out.writePixels(height);
}
