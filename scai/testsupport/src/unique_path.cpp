#include <scai/testsupport/unique_path.hpp>

#include <scai/testsupport/random_string.hpp>

#include <string>
#include <random>
#include <chrono>
#include <type_traits>

namespace
{
std::vector<std::chrono::microseconds::rep> generateSeedSequence()
{
    using std::chrono::duration_cast;
    using std::chrono::microseconds;
    using std::chrono::system_clock;
    // In the ideal case, random_device will produce truly random numbers for us, which
    // we can use to seed our random number generator. However, it is unfortunately not
    // *guaranteed* to do so, and in fact for certain implementations it is possible that
    // it will give us the same sequence of numbers every time, which would be bad for our
    // purposes. To try to alleviate this potential issue, we mix in the current time
    // to at least avoid having repeated sequences assuming that the system clock's precision
    // is sufficiently high. This inevitably leads to potential race conditions, and so it is
    // generally ill-advised, but for our purposes here it should still be sufficient,
    // as we do not need high-quality seeds.
    std::random_device device;
    const auto us_epoch = duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();

    auto r = [&device] () { return static_cast<decltype(us_epoch)>(device()); };

    return { r(), r(), r(), r(), r(), r(), r(), r(), us_epoch };
}

struct RandomGen
{
    std::mt19937 engine;
    RandomGen()
    {
        const auto seedNumbers = generateSeedSequence();
        std::seed_seq seedSequence(seedNumbers.begin(), seedNumbers.end());
        engine.seed(seedSequence);
    }
};
}

std::string scai::testsupport::unique_path(const std::string & dir, const std::string & namePrefix)
{
    static RandomGen generator;
    constexpr auto random_part_len = 8;

    const auto randomPart = randomAlphaNumericString(generator.engine, random_part_len);
    return dir + "/" + namePrefix + std::move(randomPart);
}
