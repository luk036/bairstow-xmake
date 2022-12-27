#include <doctest/doctest.h> // for ResultBuilder, CHECK, TEST_CASE

#include <bairstow/rootfinding.hpp> // for horner, initial_guess, pbairstow...
#include <utility>                  // for pair
#include <vector>                   // for vector

#include "bairstow/vector2.hpp" // for Vector2
#include "fmt/format.h"         // for print

TEST_CASE("test suppress 2") {
  auto vA = Vec2{3.0, 3.0};
  auto vA1 = Vec2{1.0, 2.0};
  const auto vri = Vec2{-2.0, 0.0};
  const auto vrj = Vec2{4.0, -5.0};

  suppress2(vA, vA1, vri, vrj);

  CHECK(vA.x() == doctest::Approx(-0.0857143));
  CHECK(vA.y() == doctest::Approx(-0.6));
}

// def test_suppress2():
//     vr0 = Vector2(-2, 0)
//     vr1 = Vector2(4, -5)
//     vr2 = Vector2(-1, 3)
//
//     vA = Vector2(3, 3)
//     vA1 = Vector2(1, 2)
//     suppress_old(vA, vA1, vr0, vr1)
//     suppress_old(vA, vA1, vr0, vr2)
//     dr_old = delta(vA, vr0, Vector2(vA1._x, -vA1._y))
//
//     vA = Vector2(3, 3)
//     vA1 = Vector2(1, 2)
//     suppress_old(vA, vA1, vr0, vr2)
//     suppress_old(vA, vA1, vr0, vr1)
//     dr_old2 = delta(vA, vr0, Vector2(vA1._x, -vA1._y))
//
//     assert dr_old2._x == approx(dr_old._x)
//     assert dr_old2._y == approx(dr_old._y)
//
//     vA = Vector2(3, 3)
//     vA1 = Vector2(1, 2)
//     suppress(vA, vA1, vr0, vr1)
//     suppress(vA, vA1, vr0, vr2)
//     dr_new = delta2(vA, vr0, vA1)
//
//     vA = Vector2(3, 3)
//     vA1 = Vector2(1, 2)
//     suppress(vA, vA1, vr0, vr2)
//     suppress(vA, vA1, vr0, vr1)
//     dr_new2 = delta2(vA, vr0, vA1)
//
//     assert dr_new._x == approx(dr_new2._x)
//     assert dr_new._y == approx(dr_new2._y)
//
//     assert dr_new._x == approx(dr_old._x)
//     assert dr_new._y == approx(dr_old._y)
