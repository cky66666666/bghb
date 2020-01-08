#include "iostream"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

double r(){
    return ((double)rand()/RAND_MAX);
}

int main(){
    for (int i = 0; i < 100; i++)
    {
        cout << ((double)rand()/RAND_MAX) << endl;
    }   
}