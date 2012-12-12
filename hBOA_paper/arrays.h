#ifndef _arrays_h_
#define _arrays_h_

class CoupleArray {

  private:

    int *x;
    int *y;
    int n;

  public:

    CoupleArray();
    ~CoupleArray();
    int add(int a, int b);
    int get(int *a, int *b, int i);
    void clear();
    int getSize();
};

#endif
