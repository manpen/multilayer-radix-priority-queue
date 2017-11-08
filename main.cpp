#include <iostream>
#include "include/multilayer_radix_pq.h"

typedef struct data
{
    int x = 0;
} DATA;

int main() {
    data x;
    multilayer_radix_pq::multilayer_radix_pq<uint64_t, DATA, 8> mlrq(4096);
    mlrq.push(pow(2,7), x);
    mlrq.push(pow(2,8), x);
    mlrq.push(pow(2,11) + pow(2,9), x);
    std::cout << "Is empty: " << mlrq.empty() << std::endl;
    return 0;
}
