
#include <iostream>

extern "C" {

double add(double a, double b) {
    return a + b;
}

void print_message(const char* message) {
    std::cout << "C++ says: " << message << std::endl;
}


}


