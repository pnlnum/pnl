#include <iostream>

extern "C" { 
    void my_hello_world ();
}

int main()
{
  my_hello_world ();
  return 0;
}
