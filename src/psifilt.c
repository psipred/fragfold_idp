#include <stdio.h>
#include <stdlib.h>

#define MIN(x,y) (((x)<(y))?(x):(y))
#define MAX(x,y) (((x)>(y))?(x):(y))

/*
Usage:

psipfilt < xxx.ss2 > xxx.out
*/

int main(int argc, char **argv)
{
   char buf[512];
   float cp, hp, ep, cutoff = 0.3;

   if (argc > 1)
       cutoff = atof(argv[1]);

   while (fgets(buf, 160, stdin))
     {
       if (sscanf(buf+8, "%f%f%f", &cp, &hp, &ep) != 3)
           continue;
       if (2*MAX(MAX(cp, hp), ep)-(cp+hp+ep)+MIN(MIN(cp, hp), ep) > cutoff)
           putchar(buf[7]);
       else if (ep > hp && cp > hp)
           putchar('e');
       else if (hp > ep && cp > ep)
           putchar('h');
       else
           putchar('?');
     }
   putchar('\n');
}
