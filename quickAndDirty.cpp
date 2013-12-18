#include <iostream>
#include <vector>
#include <climits>   //int_max

using namespace std;

//Generic sum:
template <class T> //T is a generic type
T sum(const T data[], int size, T s=0) {
    for (int i = 0; i< size; i++)
        s += data[i];
    return s;
}

//multi generic
template<class T1, class T2>
void copy(const T1 source[], T2 destination[], int size)
{
    for (int i = 0; i< size; ++i)
        destination[i] = static_cast<T2>(source[i]);    //??T1??
}

//class Point
class Point{
public:
    Point(double x=0.0, double y = 0.0): x(x), y(y) {}    //constructor
    /*
    Point(double x=0.0, double y=0.0) {
        this->x = x;
        this->y = y;
    }
    */
    double getx() {
        return this->x;
    }
    void setx(double val) {
        this->x = val;
    }
    double gety() {
        return this->y;
    }
    void sety(double val) {
        this->y = val;
    }
protected:
    double x;
    double y;
};

//list element for list
struct ListElement{
    ListElement (int n=0, ListElement* ptr=nullptr):
        d(n), next(ptr) {}
    int d;
    ListElement* next;
};

//class List:
//review this
class List{
public:
    List():head(nullptr),cursor(nullptr) {}
    void prepend(int n) {
        if (head == nullptr) {
            this->cursor = this->head = new ListElement(n,head);
        }
        else {
            head = new ListElement(n, head);
        }
    }
    int get_element() {
        return cursor->d;
    }
    void advance(){
        cursor =cursor->next;
    }
    void print() {
        ListElement* h= this->head;
        while(h != nullptr) {
            cout << h->d << ", ";
            h = h->next;
        }
        cout << "###"<< endl;
    }
protected:
  ListElement* head;
  ListElement* cursor;
};


//DIVIDE:
class Solution {
public:
public:
    int divide(long long int dividend, long long int divisor) {
        // boundary
        //if (dividend < 0)
        //    return -divide(-dividend, divisor);
        if (divisor < 0)
            return -divide(dividend, -divisor);
        if (divisor == 1)
            return dividend;
        if (dividend < divisor)
            return 0;
        if (dividend == divisor)
            return 1;
       else {
           int power = 0;
           // tocheck = divisor;
           while (divisor << power <= dividend) {
               ++power;
               //power = power << 1;
           }

           int part_result = 1 << (power - 1);
           return part_result + divide(dividend - (divisor << (power-1)), divisor);
       }
    }
};

int main()
{
//    cout << "template for sum()" << endl;
//    int a[] = {1,2,3};
//    double b[] = {2.1, 2.2, 2.3};
//    double c[] = {0,0,0};
//    //vector<int> c{3,4,1};
//    copy(a, c, 3);
//    cout << sum(a,3) << endl;
//    cout << sum(b,3) << endl;
//    cout << sum(c,3)<<endl;
    //cout << sum(c,3,4)<< endl;
    //cout << "Hello World!" <<  endl;
//    Point k(1.0);
//    cout << k.getx() << " and " << k.gety()<< endl;
//    List a, b;
//    a.prepend(9);
//    a.prepend(8);

//    a.print();

//    for(auto i=0; i< 40; ++i) {
//        b.prepend(i*i);
//    }
//    cout << "List b:" << endl;
//    b.print();
    Solution test;
    auto result = test.divide(10,3);
    cout << "max int" << INT_MAX << endl;
    cout << "result"<< result<<endl;
    return 0;
}

