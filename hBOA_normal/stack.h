#ifndef _stack_h_
#define _stack_h_

class IntStack {

 private:

  int maxSize;
  int marker;
  int size;
  int *s;

 public:

  IntStack(int max);
  ~IntStack();

  int push(int x);
  int pop();
  int setMarker();
  int goToMarker();

  int empty();
  int notEmpty();
  int full();
  int getSize();
};


#endif
