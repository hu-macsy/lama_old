#pragma once

#include <string>
#include <stdexcept>
#include <memory>
#include <utility>

namespace scai
{

namespace testsupport
{

/**
 * Represents a global temporary directory for use in tests.
 *
 * It is designed explicitly only for use in tests, with the following use pattern:
 *
 * 1. During test binary initialization, the temporary directory path is set.
 * 2. During execution of individual tests, the path can be acquired through getPath().
 *
 * Note that the class is NOT designed for concurrent or parallel access, and as such
 * it should only be used from a single thread. Note that the Boost test runner
 * is single-threaded, and so it is generally safe to use this class from tests.
 */
class GlobalTempDir
{
public:

    /**
     * Retrieves the temporary directory path.
     *
     * If this has not previously been set through a call to setPath(),
     * an std::logic_error exception is thrown.
     */
    static std::string getPath()
    {
        if (m_tempDirPath)
        {
            return *m_tempDirPath;
        }
        else
        {
            throw std::logic_error("Attempt to get the global temporary directory path "
                                   "before it has been initialized.");
        }
    }

    /**
     * Sets the temporary directory path.
     *
     * This can only be done exactly once throughout the lifetime of the program.
     * If a second attempt to set the path is made, an std::logic_error exception
     * is thrown.
     */
    static void setPath(std::string path)
    {
        if (!m_tempDirPath)
        {
            m_tempDirPath.reset(new std::string(std::move(path)));
        }
        else
        {
            throw std::logic_error("Temporary directory path has already been set. "
                                   "It must only be set exactly once in the test binary initialization.");
        }
    }

    /**
     * Sets the temporary directory path, or if the supplied string is empty,
     * set the path to a default value instead.
     *
     * The same restrictions on usage as setPath() apply to this function.
     */
    static void setPathOrDefault(std::string path)
    {
        if (path.empty())
        {
            setPath("/tmp");
        }
        else
        {
            setPath(path);
        }
    }

private:
    // Using unique_ptr as a poor man's optional...
    static std::unique_ptr<std::string> m_tempDirPath;
};

} // namespace testsupport

} // namespace scai
