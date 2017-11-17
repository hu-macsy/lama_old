#pragma once

#include <string>
#include <scai/hmemo/Context.hpp>

namespace scai
{

namespace testsupport
{

namespace detail
{

/**
 * Find a suitable name for a test suite which describes
 * the context used. Note that the naming is chosen such that
 * related environments are collected when sorted.
 */
std::string adaptTestSuiteNameToEnv(const std::string & name, const scai::hmemo::Context & context)
{
    std::string prefix;
    switch (context.getType())
    {
        case scai::common::context::Host:
            prefix = "~Host ";
            break;
        case scai::common::context::CUDA:
            prefix = "~CUDA ";
            break;
        default:
            std::cerr << "Unsupported context type. Can not create appropriate test suite name." << std::endl;
            throw std::runtime_error("Unsupported context type.");
    }

    return prefix + name;
}

std::string suiteNameForFile(const std::string & name, const scai::hmemo::Context & context)
{
    // TODO: Context
    std::stringstream filename;
    filename << name;

    switch (context.getType())
    {
        case scai::common::context::Host:
            filename << "_host";
            break;
        case scai::common::context::CUDA:
            filename << "_cuda";
            break;
        default:
            std::cerr << "Unsupported context type. Can not create appropriate test filename." << std::endl;
            throw std::runtime_error("Unsupported context type.");
    }

    return filename.str();
}

} // namespace detail

} // namespace testsupport

} // namespace scai
