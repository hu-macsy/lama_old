#include <memory>
#include <iostream>
 
struct Foo
{
    Foo()  { std::cout << "Foo...\n";  }
    ~Foo() { std::cout << "~Foo...\n"; }
};
 
struct D
{
    void operator()(Foo* p) const
    {
        delete p;
    }
};
 
int main()
{
    {
        // constructor with no managed object
        std::shared_ptr<Foo> sh1;
    }
 
    {
        // constructor with object
        std::shared_ptr<Foo> sh2(new Foo);
        std::shared_ptr<Foo> sh3(sh2);
        std::cout << sh2.use_count() << '\n';
        std::cout << sh3.use_count() << '\n';
    }
 
    {
        // constructor with object and deleter
        std::shared_ptr<Foo> sh4(new Foo, D());
    }

    return 0;
}