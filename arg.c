#include <stdio.h>
typedef struct test
{
  unsigned char  *** t_arr; //Actual data array
  union {
    int a;
    float b;
    double c;
  };
} test;

void main(void) {
  test jo;
  jo.a  =5;
  jo.b = 3.42323;
  jo.c = 3.12121;
  jo.a = 100;
  printf("%d", jo.a);

}
