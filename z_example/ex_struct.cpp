#include <iostream>
using namespace std;

struct MyStruct {
  int x;
  int y;
};

class MyClass {
  public:
    MyStruct myArray[10];
    void printMyArray() {
      for (int i = 0; i < 10; i++) {
        cout << "x: " << myArray[i].x << ", y: " << myArray[i].y << endl;
      }
    }
};

int main() {
  MyClass myClass;
  for (int i = 0; i < 10; i++) {
    myClass.myArray[i].x = i;
    myClass.myArray[i].y = i * 2;
  }
  myClass.printMyArray();
}