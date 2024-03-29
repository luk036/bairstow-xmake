#include <doctest/doctest.h> // for ResultBuilder, CHECK, Expr...
// #include <__config>                        // for std
#include <bairstow/robin.hpp> // for Robin, Robin<>::iterable_w...
#include <cinttypes>          // for uint8_t
#include <utility>            // for pair

using namespace std;

TEST_CASE("Test Robin") {
  fun::Robin<uint8_t> rr(6U);
  auto count = 0U;
  for (auto _i : rr.exclude(2)) {
    static_assert(sizeof _i >= 0, "make compiler happy");
    count += 1;
  }
  CHECK(count == 5);
}
