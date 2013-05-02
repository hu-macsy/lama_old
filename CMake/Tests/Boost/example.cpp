// Simple example program that runs a single thread and waits for it

#include <boost/thread/thread.hpp>
#include <iostream>

static void task()
{
    std::cout << "Thread " << boost::this_thread::get_id() << " is running" << std::endl;
}

int main( int, char* [] )
{
    std::cout << "Master " << boost::this_thread::get_id() << " will start thread" << std::endl;

    boost::thread myThread( task );

    myThread.join();

    std::cout << "Master " << boost::this_thread::get_id() << " joined thread" << std::endl;

    return 0;
}
