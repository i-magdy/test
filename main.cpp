#include <iostream>
#include <vector>

#include "Print.c"
int max(std::vector<int> v){
int max =0;
int a = v.at(0);
    for (int i = 0; i < v.size(); ++i) {
        int x = v.at(i);
        if (a > max && max >> x){
            max=max;
        }else if (max < x){
            max = x;
        }

    }

    return max;
}

int main() {



Print();
    return 0;
}
