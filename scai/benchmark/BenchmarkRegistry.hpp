/**
 * @file BenchmarkRegistry.h
 *
 * @license
 * Copyright (c) 2011
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief BenchmarkRegistry.hpp
 * @author jiri
 * @date 06.04.2011
 * $Id$
 */
#include <memory>
#include <string>
#include <map>
#include <iostream>

#include <scai/common/config.hpp>
#include <scai/benchmark/Benchmark.hpp>
#include <scai/benchmark/macros.hpp>

/**
 * @brief The namespace bf holds everything of the benchmark framework
 */
namespace bf
{
class COMMON_DLL_IMPORTEXPORT AbstractBenchmarkCreator
{
public:
    /**
     * @brief Destructor destroys this object and frees all inner ressources, if
     * necessary.
     */
    virtual ~AbstractBenchmarkCreator();
    /**
     * @brief The function to create the benchmark.
     */
    virtual Benchmark* create() const =0;
    /**
     * @brief The function to create the benchmark with arguments.
     * @param[in] args The arguments passed to the benchmark.
     */
    virtual Benchmark* create( const std::string& args ) const =0;
};

class COMMON_DLL_IMPORTEXPORT BenchmarkRegistry
{
public:
    typedef std::map<std::string,AbstractBenchmarkCreator*>::const_iterator const_iterator;

    /**
     * @brief Destructor destroys this object and frees all inner ressources, if
     *        necessary.
     */
    virtual ~BenchmarkRegistry();

    /**
     * @brief Returns the registry. If not yet existing, it will be created.
     * @return The registry.
     */
    static BenchmarkRegistry& getRegistry();

    /**
     * @brief Releases all resource of the registry. Once this method is called it is no longer possible to
     *        restore the orignial state.
     */
    static void freeRegistry();

    /**
     * @brief Creates benchmark with given id.
     * @param[in] id The id of the benchmark to be created.
     * @return The benchmark in an autopointer.
     * @throws Exception, if the id could not be resolved to a benchmark.
     *
     * @todo should return no auto_ptr a managment object that calls destroyBenchmark instead of delete is necessary (TODO)
     */
    std::auto_ptr<Benchmark>
    createBenchmark( const std::string& id ) const;

    /**
     * @brief Creates benchmark with given id.
     * @param[in] id   The id of the benchmark to be created.
     * @param[in] args The arguments, passed to the benchmark.
     * @return The benchmark in an autopointer.
     * @throws Exception, if the id could not be resolved to a benchmark.
     *
     * @todo should return no auto_ptr a managment object that calls destroyBenchmark instead of delete is necessary (TODO)
     */
    std::auto_ptr<Benchmark>
    createBenchmark( const std::string& id, const std::string& args ) const;

    void destroyBenchmark( Benchmark* bench ) const;

    /**
     * @brief Adds a benchmark to the registry.
     * @param[in] id        The unique id of the benchmark.
     * @param[in] creator   The function to create the benchmark.
     */
    void add( const std::string& id, AbstractBenchmarkCreator* creator );

    bool has( const std::string& id ) const;

    /**
     * @brief Returns iterator to beginning.
     * @return Iterator to beginning.
     */
    const_iterator begin() const;

    /**
     * @brief Returns iterator to end.
     * @return Iterator to end.
     */
    const_iterator end() const;

private:
    /**
     * @brief Returns the pointer to an AbstractBenchmarkCreator to a given ID.
     * @param[in] id The ID of the benchmark, created by the
     *               AbstractBenchmarkCreator.
     * @throws Exception, if the id could not be resolved to a benchmark.
     */
    AbstractBenchmarkCreator* getCreatorById( const std::string& id ) const;
    /**
     * @brief Default constructor.
     */
    BenchmarkRegistry();
    /**
     * @brief Copy constructor.
     */
    BenchmarkRegistry( const BenchmarkRegistry& other );
    static BenchmarkRegistry* mInstance;
    std::map<std::string,AbstractBenchmarkCreator*> mCreatorMap;

    static bool isFreed;

    class CGuard
    {
    public:
        /**
         * @brief Default constructor.
         */
        CGuard();
        /**
         * @brief Destructor destroys array of BenchmarkRegistry, if initialized.
         */
        ~CGuard();
    };
    friend class CGuard;
};

template<typename T>
class BenchmarkCreator: public AbstractBenchmarkCreator
{
public:
    virtual T* create() const
    {
        return new T();
    }

    virtual T* create( const std::string& args ) const
    {
        return new T( args );
    }
};

template<typename T>
class BenchmarkRegistration
{
public:
    BenchmarkRegistration( const std::string& id )
    {
        BenchmarkRegistry& reg = BenchmarkRegistry::getRegistry();
        reg.add( id, new BenchmarkCreator<T>() );
    }
};

/**
 * @brief Register the type.
 */
#define LAMA_BENCHMARK_REGISTRATION(type)                                            \
    static bf::BenchmarkRegistration<type >                                  \
    LAMA_UNIQUE_NAME( benchRegObj, post )(type::id())
/**
 * @brief Register the type.
 */
#define LAMA_BENCHMARK_REGISTRATION2(type,post)                                      \
    static bf::BenchmarkRegistration<type >                                  \
    LAMA_UNIQUE_NAME( benchRegObj, post )(type::id())

}

