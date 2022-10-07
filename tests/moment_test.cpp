#include <gtest/gtest.h>
#include <moments.h>
#include <moments_image.h>

#include <vector>
#include "moment_test_data.h"

struct Spectrum {
    // const std::vector<double>& wavelengths;
    bool isReflective;
    double wavelength_start, wavelength_end;
    int n_wavelengths;
    const std::vector<double>& values;
    const char* name;
};


class MomentTest: public ::testing::TestWithParam<Spectrum> {
protected:
    void SetUp() override
    {
        const Spectrum& s = GetParam();
        isReflective = s.isReflective;
        signal = s.values;
        linspace(s.wavelength_start, s.wavelength_end, s.n_wavelengths, wavelengths);
        wavelengths_to_phases(wavelengths, phases);
        compute_moments(phases, s.values, phases.size(), moments);
    }

    bool isReflective;
    std::vector<double> signal;
    std::vector<double> wavelengths;
    std::vector<double> phases;
    std::vector<double> moments;
};


TEST_P(MomentTest, Sizes)
{
    ASSERT_EQ(wavelengths.size(), phases.size());
    ASSERT_EQ(wavelengths.size(), moments.size());
}


TEST_P(MomentTest, BackDensity)
{
    std::vector<double> signal_back;

    if (!isReflective) {
        bounded_compute_density_lagrange(phases, moments, signal_back);
    } else {
        compute_density(phases, moments, signal_back);
    }
    
    ASSERT_EQ(phases.size(), signal_back.size());
    ASSERT_EQ(signal.size(), signal_back.size());

    for (size_t i = 0; i < signal.size(); i++) {
        EXPECT_NEAR(signal[i], signal_back[i], 1e-1);
    }
}


TEST_P(MomentTest, BackPointBased)
{
    std::vector<double> signal_back;
    compute_density_image(phases, moments, 1, 1, moments.size(), signal_back);

    ASSERT_EQ(signal.size(), signal_back.size());

    for (size_t i = 0; i < signal.size(); i++) {
        EXPECT_NEAR(signal[i], signal_back[i], 1e-10);
    }
}


TEST_P(MomentTest, BoundedCompress)
{
    if (!isReflective) {
        // TODO properly implement the test for reflective spectra
        GTEST_SKIP();
    }

    std::vector<double> compressed_moments;
    std::vector<double> inflated_moments;

    bounded_compress_moments(moments, compressed_moments);
    ASSERT_EQ(moments.size(), compressed_moments.size());

    bounded_decompress_moments(compressed_moments, inflated_moments);
    ASSERT_EQ(compressed_moments.size(), inflated_moments.size());
    
    for (size_t i = 0; i < moments.size(); i++) {
        EXPECT_NEAR(moments[i], inflated_moments[i], 1e-8);
    }
}


TEST_P(MomentTest, UnboundedCompress)
{
    std::vector<double> compressed_moments;
    std::vector<double> inflated_moments;

    unbounded_compress_moments(moments, compressed_moments);
    ASSERT_EQ(moments.size(), compressed_moments.size());

    unbounded_decompress_moments(compressed_moments, inflated_moments);
    ASSERT_EQ(compressed_moments.size(), inflated_moments.size());
    
    for (size_t i = 0; i < moments.size(); i++) {
        EXPECT_NEAR(moments[i], inflated_moments[i], 1e-8);
    }
}


TEST_P(MomentTest, UnboundedToBoundedCompress)
{
    std::vector<double> compressed_moments;
    std::vector<double> inflated_moments;

    unbounded_to_bounded_compress_moments(moments, compressed_moments);
    ASSERT_EQ(moments.size(), compressed_moments.size());

    unbounded_to_bounded_decompress_moments(compressed_moments, inflated_moments);
    ASSERT_EQ(compressed_moments.size(), inflated_moments.size());
    
    for (size_t i = 0; i < moments.size(); i++) {
        EXPECT_NEAR(moments[i], inflated_moments[i], 1e-8);
    }
}


Spectrum spectrum_reflectance = {true, reflectance_start, reflectance_end, reflectance_n_samples, reflectance, "Reflectance"};
Spectrum spectrum_emission    = {false, emission_start, emission_end, emission_n_samples, emission, "Emission"};

// Macbeth
Spectrum s_macbeth_dark_skin     = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_dark_skin, "DarkSkin"};
Spectrum s_macbeth_light_skin    = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_light_skin, "LightSkin"};
Spectrum s_macbeth_blue_sky      = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_blue_sky, "BlueSky"};
Spectrum s_macbeth_foliage       = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_foliage, "Foliage"};
Spectrum s_macbeth_blue_flower   = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_blue_flower, "BlueFlower"};
Spectrum s_macbeth_bluish_green  = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_bluish_green, "BluishGreen"};

Spectrum s_macbeth_orange        = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_orange, "Orange"};
Spectrum s_macbeth_purplish_blue = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_purplish_blue, "PurpulishBlue"};
Spectrum s_macbeth_moderate_red  = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_moderate_red, "ModerateRed"};
Spectrum s_macbeth_purple        = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_purple, "Purple"};
Spectrum s_macbeth_yellow_green  = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_yellow_green, "YellowGreen"};
Spectrum s_macbeth_orange_yellow = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_orange_yellow, "OrangeYellow"};

Spectrum s_macbeth_blue          = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_blue, "Blue"};
Spectrum s_macbeth_green         = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_green, "Green"};
Spectrum s_macbeth_red           = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_red, "Red"};
Spectrum s_macbeth_yellow        = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_yellow, "Yellow"};
Spectrum s_macbeth_magenta       = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_magenta, "Magenta"};
Spectrum s_macbeth_cyan          = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_cyan, "Cyan"};

Spectrum s_macbeth_white         = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_white, "White"};
Spectrum s_macbeth_neutral_8     = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_neutral_8, "Neutral8"};
Spectrum s_macbeth_neutral_6_5   = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_neutral_6_5, "Neutral65"};
Spectrum s_macbeth_neutral_5     = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_neutral_5, "Neutral5"};
Spectrum s_macbeth_neutral_3_5   = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_neutral_3_5, "Neutral35"};
Spectrum s_macbeth_neutral_black = {true, macbeth_start, macbeth_end, macbeth_n_samples, macbeth_neutral_black, "Black"};

// Standard illuminants
Spectrum s_FL3_13                               = { false, FL3_13_start, FL3_13_end, FL3_13_n_samples, FL3_13, "FL3_13" };
Spectrum s_FL3_11                               = { false, FL3_11_start, FL3_11_end, FL3_11_n_samples, FL3_11, "FL3_11" };
Spectrum s_ISO_7589_Sensitometric_Studio_Tungsten = { false, ISO_7589_Sensitometric_Studio_Tungsten_start, ISO_7589_Sensitometric_Studio_Tungsten_end, ISO_7589_Sensitometric_Studio_Tungsten_n_samples, ISO_7589_Sensitometric_Studio_Tungsten, "ISO_7589_Sensitometric_Studio_Tungsten" };
Spectrum s_LED_V1                               = { false, LED_V1_start, LED_V1_end, LED_V1_n_samples, LED_V1, "LED_V1" };
Spectrum s_FL8                                  = { false, FL8_start, FL8_end, FL8_n_samples, FL8, "FL8" };
Spectrum s_ISO_7589_Studio_Tungsten             = { false, ISO_7589_Studio_Tungsten_start, ISO_7589_Studio_Tungsten_end, ISO_7589_Studio_Tungsten_n_samples, ISO_7589_Studio_Tungsten, "ISO_7589_Studio_Tungsten" };
Spectrum s_ID50                                 = { false, ID50_start, ID50_end, ID50_n_samples, ID50, "ID50" };
Spectrum s_FL11                                 = { false, FL11_start, FL11_end, FL11_n_samples, FL11, "FL11" };
Spectrum s_FL2                                  = { false, FL2_start, FL2_end, FL2_n_samples, FL2, "FL2" };
Spectrum s_FL3_3                                = { false, FL3_3_start, FL3_3_end, FL3_3_n_samples, FL3_3, "FL3_3" };
Spectrum s_FL5                                  = { false, FL5_start, FL5_end, FL5_n_samples, FL5, "FL5" };
Spectrum s_LED_V2                               = { false, LED_V2_start, LED_V2_end, LED_V2_n_samples, LED_V2, "LED_V2" };
Spectrum s_FL3_1                                = { false, FL3_1_start, FL3_1_end, FL3_1_n_samples, FL3_1, "FL3_1" };
Spectrum s_ISO_7589_Sensitometric_Printer       = { false, ISO_7589_Sensitometric_Printer_start, ISO_7589_Sensitometric_Printer_end, ISO_7589_Sensitometric_Printer_n_samples, ISO_7589_Sensitometric_Printer, "ISO_7589_Sensitometric_Printer" };
Spectrum s_LED_B2                               = { false, LED_B2_start, LED_B2_end, LED_B2_n_samples, LED_B2, "LED_B2" };
Spectrum s_FL4                                  = { false, FL4_start, FL4_end, FL4_n_samples, FL4, "FL4" };
Spectrum s_C                                    = { false, C_start, C_end, C_n_samples, C, "C" };
Spectrum s_E                                    = { false, E_start, E_end, E_n_samples, E, "E" };
Spectrum s_FL10                                 = { false, FL10_start, FL10_end, FL10_n_samples, FL10, "FL10" };
Spectrum s_FL3_15                               = { false, FL3_15_start, FL3_15_end, FL3_15_n_samples, FL3_15, "FL3_15" };
Spectrum s_FL3_6                                = { false, FL3_6_start, FL3_6_end, FL3_6_n_samples, FL3_6, "FL3_6" };
Spectrum s_FL3_5                                = { false, FL3_5_start, FL3_5_end, FL3_5_n_samples, FL3_5, "FL3_5" };
Spectrum s_FL12                                 = { false, FL12_start, FL12_end, FL12_n_samples, FL12, "FL12" };
Spectrum s_D60                                  = { false, D60_start, D60_end, D60_n_samples, D60, "D60" };
Spectrum s_D75                                  = { false, D75_start, D75_end, D75_n_samples, D75, "D75" };
Spectrum s_D55                                  = { false, D55_start, D55_end, D55_n_samples, D55, "D55" };
Spectrum s_FL1                                  = { false, FL1_start, FL1_end, FL1_n_samples, FL1, "FL1" };
Spectrum s_HP4                                  = { false, HP4_start, HP4_end, HP4_n_samples, HP4, "HP4" };
Spectrum s_HP1                                  = { false, HP1_start, HP1_end, HP1_n_samples, HP1, "HP1" };
Spectrum s_FL3                                  = { false, FL3_start, FL3_end, FL3_n_samples, FL3, "FL3" };
Spectrum s_FL3_4                                = { false, FL3_4_start, FL3_4_end, FL3_4_n_samples, FL3_4, "FL3_4" };
Spectrum s_LED_BH1                              = { false, LED_BH1_start, LED_BH1_end, LED_BH1_n_samples, LED_BH1, "LED_BH1" };
Spectrum s_LED_RGB1                             = { false, LED_RGB1_start, LED_RGB1_end, LED_RGB1_n_samples, LED_RGB1, "LED_RGB1" };
Spectrum s_FL3_8                                = { false, FL3_8_start, FL3_8_end, FL3_8_n_samples, FL3_8, "FL3_8" };
Spectrum s_HP5                                  = { false, HP5_start, HP5_end, HP5_n_samples, HP5, "HP5" };
Spectrum s_A                                    = { false, A_start, A_end, A_n_samples, A, "A" };
Spectrum s_D65                                  = { false, D65_start, D65_end, D65_n_samples, D65, "D65" };
Spectrum s_FL6                                  = { false, FL6_start, FL6_end, FL6_n_samples, FL6, "FL6" };
Spectrum s_HP3                                  = { false, HP3_start, HP3_end, HP3_n_samples, HP3, "HP3" };
Spectrum s_ISO_7589_Photographic_Daylight       = { false, ISO_7589_Photographic_Daylight_start, ISO_7589_Photographic_Daylight_end, ISO_7589_Photographic_Daylight_n_samples, ISO_7589_Photographic_Daylight, "ISO_7589_Photographic_Daylight" };
Spectrum s_LED_B5                               = { false, LED_B5_start, LED_B5_end, LED_B5_n_samples, LED_B5, "LED_B5" };
Spectrum s_HP2                                  = { false, HP2_start, HP2_end, HP2_n_samples, HP2, "HP2" };
Spectrum s_FL3_10                               = { false, FL3_10_start, FL3_10_end, FL3_10_n_samples, FL3_10, "FL3_10" };
Spectrum s_B                                    = { false, B_start, B_end, B_n_samples, B, "B" };
Spectrum s_FL3_14                               = { false, FL3_14_start, FL3_14_end, FL3_14_n_samples, FL3_14, "FL3_14" };
Spectrum s_FL7                                  = { false, FL7_start, FL7_end, FL7_n_samples, FL7, "FL7" };
Spectrum s_ISO_7589_Sensitometric_Daylight      = { false, ISO_7589_Sensitometric_Daylight_start, ISO_7589_Sensitometric_Daylight_end, ISO_7589_Sensitometric_Daylight_n_samples, ISO_7589_Sensitometric_Daylight, "ISO_7589_Sensitometric_Daylight" };
Spectrum s_FL3_7                                = { false, FL3_7_start, FL3_7_end, FL3_7_n_samples, FL3_7, "FL3_7" };
Spectrum s_FL3_12                               = { false, FL3_12_start, FL3_12_end, FL3_12_n_samples, FL3_12, "FL3_12" };
Spectrum s_ID65                                 = { false, ID65_start, ID65_end, ID65_n_samples, ID65, "ID65" };
Spectrum s_D50                                  = { false, D50_start, D50_end, D50_n_samples, D50, "D50" };
Spectrum s_FL3_2                                = { false, FL3_2_start, FL3_2_end, FL3_2_n_samples, FL3_2, "FL3_2" };
Spectrum s_ISO_7589_Photoflood                  = { false, ISO_7589_Photoflood_start, ISO_7589_Photoflood_end, ISO_7589_Photoflood_n_samples, ISO_7589_Photoflood, "ISO_7589_Photoflood" };
Spectrum s_FL9                                  = { false, FL9_start, FL9_end, FL9_n_samples, FL9, "FL9" };
Spectrum s_LED_B1                               = { false, LED_B1_start, LED_B1_end, LED_B1_n_samples, LED_B1, "LED_B1" };
Spectrum s_ISO_7589_Sensitometric_Photoflood    = { false, ISO_7589_Sensitometric_Photoflood_start, ISO_7589_Sensitometric_Photoflood_end, ISO_7589_Sensitometric_Photoflood_n_samples, ISO_7589_Sensitometric_Photoflood, "ISO_7589_Sensitometric_Photoflood" };
Spectrum s_LED_B3                               = { false, LED_B3_start, LED_B3_end, LED_B3_n_samples, LED_B3, "LED_B3" };
Spectrum s_FL3_9                                = { false, FL3_9_start, FL3_9_end, FL3_9_n_samples, FL3_9, "FL3_9" };
Spectrum s_LED_B4                               = { false, LED_B4_start, LED_B4_end, LED_B4_n_samples, LED_B4, "LED_B4" };


INSTANTIATE_TEST_SUITE_P(
    Reflectance, MomentTest,
    ::testing::Values(
        spectrum_reflectance,
        // spectrum_emission,
        // Macbeth colour checker
        s_macbeth_dark_skin,
        s_macbeth_light_skin,
        s_macbeth_blue_sky,
        s_macbeth_foliage,
        s_macbeth_blue_flower,
        s_macbeth_bluish_green,
        //
        s_macbeth_orange,
        s_macbeth_purplish_blue,
        s_macbeth_moderate_red,
        s_macbeth_purple,
        s_macbeth_yellow_green,
        s_macbeth_orange_yellow,
        //
        s_macbeth_blue,
        s_macbeth_green,
        s_macbeth_red,
        s_macbeth_yellow,
        s_macbeth_magenta,
        s_macbeth_cyan,
        //
        s_macbeth_white,
        s_macbeth_neutral_8,
        s_macbeth_neutral_6_5,
        s_macbeth_neutral_5,
        s_macbeth_neutral_3_5,
        s_macbeth_neutral_black,
        // Standard Illuminants
        s_FL3_13,
        s_FL3_11,
        s_ISO_7589_Sensitometric_Studio_Tungsten,
        s_LED_V1,
        s_FL8,
        s_ISO_7589_Studio_Tungsten,
        s_ID50,
        s_FL11,
        s_FL2,
        s_FL3_3,
        s_FL5,
        s_LED_V2,
        s_FL3_1,
        s_ISO_7589_Sensitometric_Printer,
        s_LED_B2,
        s_FL4,
        s_C,
        s_E,
        s_FL10,
        s_FL3_15,
        s_FL3_6,
        s_FL3_5,
        s_FL12,
        s_D60,
        s_D75,
        s_D55,
        s_FL1,
        s_HP4,
        s_HP1,
        s_FL3,
        s_FL3_4,
        s_LED_BH1,
        s_LED_RGB1,
        s_FL3_8,
        s_HP5,
        s_A,
        s_D65,
        s_FL6,
        s_HP3,
        s_ISO_7589_Photographic_Daylight,
        s_LED_B5,
        s_HP2,
        s_FL3_10,
        s_B,
        s_FL3_14,
        s_FL7,
        s_ISO_7589_Sensitometric_Daylight,
        s_FL3_7,
        s_FL3_12,
        s_ID65,
        s_D50,
        s_FL3_2,
        s_ISO_7589_Photoflood,
        s_FL9,
        s_LED_B1,
        s_ISO_7589_Sensitometric_Photoflood,
        s_LED_B3,
        s_FL3_9,
        s_LED_B4
        ),
    [](const ::testing::TestParamInfo<Spectrum> &info)
    {
        return info.param.name;
    });
