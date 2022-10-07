#include <gtest/gtest.h>
#include <Util.h>

#include "util_test_data.h"

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


TEST(Util, Interp)
{
    ASSERT_FLOAT_EQ(Util::interp(0., -3., 3., -10., 0.), -5.0);

    ASSERT_FLOAT_EQ(Util::interp(3., 3., 10., 5., 10.), 5.0);    
    ASSERT_FLOAT_EQ(Util::interp(10., 3., 10., 5., 10.), 10.0);    
}