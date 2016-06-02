#include <iostream>
#include <memory>
 
void observe( std::weak_ptr<int> weak ) 
{
    if ( std::shared_ptr<int> observe = weak.lock() )
    {
        std::cout << "\tobserve() able to lock weak_ptr<>, value=" << *observe << "\n";
    }
    else
    {
        std::cout << "\tobserve() unable to lock weak_ptr<>\n";
    }
}
 
int main()
{
    std::weak_ptr<int> weak;
    std::cout << "weak_ptr<> not yet initialized\n";
    observe( weak );
 
    {
        std::weak_ptr<int> shared = std::make_shared<int>( 42 );
        weak = shared;
        std::cout << "weak_ptr<> initialized with shared_ptr.\n";
        observe( weak );
    }
 
    std::cout << "shared_ptr<> has been destructed due to scope exit.\n";
    observe( weak );

    return 0;
}