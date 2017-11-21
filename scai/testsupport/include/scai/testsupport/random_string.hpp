#pragma once

#include <random>
#include <string>
#include <cassert>

namespace scai
{

namespace testsupport
{

template <typename RandomEngine>
std::string randomAlphaNumericString(RandomEngine & engine, size_t len)
{
    const static std::string allowed_chars =
        "01234567890"
        "abcdefghjiklmnopqrstuvwxyz"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

    // Need to use unsigned int instead of size_t because it is not guaranteed
    // that a distribution for size_t exists
    std::uniform_int_distribution<unsigned int> distribution(
        static_cast<unsigned int>(0),
        static_cast<unsigned int>(allowed_chars.size() - 1)
    );

    std::string s;
    s.reserve(len);
    for (size_t i = 0; i < len; ++i)
    {
        const auto random_index = distribution(engine);
        assert(random_index < allowed_chars.size());
        const auto random_char = allowed_chars[static_cast<size_t>(random_index)];
        s.push_back(random_char);
    }

    return s;
}

} // namespace testsupport

} // namespace scai
