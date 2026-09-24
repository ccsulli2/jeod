/*
 * gravity_controls_ut.cc
 */

#include "environment/gravity/include/gravity_controls.hh"
#include "environment/gravity/include/gravity_integ_frame.hh"
#include "environment/gravity/include/gravity_manager.hh"
#include "environment/gravity/include/gravity_source.hh"
#include "environment/ephemerides/ephem_interface/include/ephem_ref_frame.hh"
#include "message_handler_mock.hh"

#include "gmock/gmock.h"
#include <cmath>
#include "gtest/gtest.h"
using testing::_;
using testing::AnyNumber;
using testing::Mock;

using namespace jeod;

TEST(GravityControls, initialize_control) {}

TEST(GravityControls, reset_control) {}

TEST(GravityControls, gravitation) {}

class TestGravityControls : public GravityControls
{
public:
    using GravityControls::calc_spherical;

    void set_manager(GravityManager & manager)
    {
        grav_manager = &manager;
    }

protected:
    void calc_nonspherical(const double[3],
                           const double[3],
                           const GravityIntegFrame &,
                           double[3],
                           double[3][3],
                           double &) override
    {
    }
};

TEST(GravityControls, embary_third_body)
{
    MockMessageHandler mockMessageHandler;
    EXPECT_CALL(mockMessageHandler, process_message(_, _, _, _, _, _, _)).Times(AnyNumber());

    EphemerisRefFrame ssbary;
    EphemerisRefFrame sun_frame;
    EphemerisRefFrame embary;
    EphemerisRefFrame earth_frame;
    EphemerisRefFrame moon_frame;
    embary.set_name("EMBary.inertial");
    ssbary.add_child(sun_frame);
    ssbary.add_child(embary);
    embary.add_child(earth_frame);
    embary.add_child(moon_frame);
    embary.state.trans.position[0] = 10.0;
    earth_frame.state.trans.position[0] = -1.0;
    moon_frame.state.trans.position[0] = 3.0;

    GravitySource sun;
    GravitySource earth;
    GravitySource moon;
    sun.name = "Sun";
    sun.mu = 100.0;
    sun.inertial = &sun_frame;
    earth.name = "Earth";
    earth.mu = 3.0;
    earth.inertial = &earth_frame;
    moon.name = "Moon";
    moon.mu = 1.0;
    moon.inertial = &moon_frame;
    GravityManager manager;
    manager.add_grav_source(sun);
    manager.add_grav_source(earth);
    manager.add_grav_source(moon);

    TestGravityControls control;
    control.body = &sun;
    control.set_manager(manager);

    GravityIntegFrame frame;
    frame.ref_frame = &embary;
    frame.is_third_body = true;
    frame.pos[0] = 10.0; // Sun to EMBary

    const double integ_pos[3] = {2.0, 1.0, 0.0};
    const double posn[3] = {12.0, 1.0, 0.0}; // Sun to vehicle
    const double direct_x = -100.0 * 12.0 / std::pow(145.0, 1.5);
    const double indirect_x = 100.0 * (0.75 / 81.0 + 0.25 / 169.0);

    for(bool battin : {false, true})
    {
        control.battin_method = battin;
        double accel[3] = {};
        double gradient[3][3] = {};
        double potential = 0.0;
        control.calc_spherical(integ_pos, posn, frame, accel, gradient, potential);
        EXPECT_NEAR(accel[0], direct_x + indirect_x, 1e-12);
        EXPECT_NEAR(accel[1], -100.0 / std::pow(145.0, 1.5), 1e-12);
        EXPECT_DOUBLE_EQ(accel[2], 0.0);
        EXPECT_NEAR(potential, 100.0 / std::sqrt(145.0), 1e-12);
    }

    // An ordinary Earth-centered frame still uses the point-origin offset.
    frame.ref_frame = &earth_frame;
    frame.pos[0] = 9.0;
    double accel[3] = {};
    double gradient[3][3] = {};
    double potential = 0.0;
    control.battin_method = false;
    control.calc_spherical(integ_pos, posn, frame, accel, gradient, potential);
    EXPECT_NEAR(accel[0], direct_x + 100.0 / 81.0, 1e-12);
}

TEST(GravityControls, calc_relativistic) {}
