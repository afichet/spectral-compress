#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <Util.h>

#include "util_test_data.h"


TEST(Util, Interp)
{
    ASSERT_FLOAT_EQ(Util::interp(0., -3., 3., -10., 0.), -5.0);

    ASSERT_FLOAT_EQ(Util::interp(3., 3., 10., 5., 10.), 5.0);    
    ASSERT_FLOAT_EQ(Util::interp(10., 3., 10., 5., 10.), 10.0);    
}


TEST(Util, DeltaE2000)
{
    for (int i = 0; i < 33; i++) {
        double Lab1[3], Lab2[3];
        memcpy(Lab1, &dE_ref[7 * i + 0], 3 * sizeof(double));
        memcpy(Lab2, &dE_ref[7 * i + 3], 3 * sizeof(double));
        
        const double ref = dE_ref[7 * i + 6];

        const double v1 = Util::deltaE2000(Lab1, Lab2);
        const double v2 = Util::deltaE2000(Lab2, Lab1);

        EXPECT_NEAR(ref, v1, 1e-4);
        EXPECT_NEAR(ref, v2, 1e-4);
    }
}


TEST(Util, Linespace)
{
    const std::vector<int> ref_idx = {1, 2, 3, 4, 5, 6, 7};
    
    std::vector<int> idx;

    Util::linspace(1, 7, ref_idx.size(), idx);

    ASSERT_EQ(ref_idx.size(), idx.size());

    for (size_t i = 0; i < idx.size(); i++) {
        EXPECT_EQ(ref_idx[i], idx[i]);
    }
}


TEST(Util, SplitExntension)
{
    
    std::string base, ext;

    Util::split_extension(".test", base, ext);
    EXPECT_THAT(base, ::testing::StrEq(""));
    EXPECT_THAT(ext , ::testing::StrEq(".test"));

    Util::split_extension("", base, ext);
    EXPECT_THAT(base, ::testing::StrEq(""));
    EXPECT_THAT(ext , ::testing::StrEq(""));

    Util::split_extension("test.", base, ext);
    EXPECT_THAT(base, ::testing::StrEq("test"));
    EXPECT_THAT(ext , ::testing::StrEq("."));

    Util::split_extension("test.jxl", base, ext);
    EXPECT_THAT(base, ::testing::StrEq("test"));
    EXPECT_THAT(ext , ::testing::StrEq(".jxl"));
}