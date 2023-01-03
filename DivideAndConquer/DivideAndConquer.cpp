#include <iostream>
#include <stack> 
#include <stdlib.h> 
#include <omp.h>
#include <time.h>
#include <chrono>
#include<vector>
#include <set>
#include <algorithm>
#include <execution>
#include "MonotoneChainCopy.cpp"
using namespace std;
using namespace std::chrono;

// A divide and conquer program to find convex
// hull of a given set of points.

//#define llu long long int
#define SIZE 5000000
/*

};*/

class Timer {
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1, 1000000>> second_;
    std::chrono::time_point<clock_> beg_;
    const char* header;
public:
    Timer(const char* header = "") : beg_(clock_::now()), header(header) {}
    ~Timer() {
        double e = elapsed();
        cout << header << ": " << e / 1000000 << " seconds" << endl;
    }
    void reset() {
        beg_ = clock_::now();
    }
    double elapsed() const {
        return std::chrono::duration_cast<second_>(clock_::now() - beg_).count();
    }
};



// stores the centre of polygon (It is made
// global because it is used in compare function)
Point mid;
// determines the quadrant of a point
// (used in compare())
int quad(Point p)
{
    if (p.x >= 0 && p.y >= 0)
        return 1;
    if (p.x <= 0 && p.y >= 0)
        return 2;
    if (p.x <= 0 && p.y <= 0)
        return 3;
    return 4;
}
// Checks whether the line is crossing the polygon
int orientation(Point a, Point b,
    Point c)
{
    int res = (b.y - a.y) * (c.x - b.x) -
        (c.y - b.y) * (b.x - a.x);
    if (res == 0)
        return 0;
    if (res > 0)
        return 1;
    return -1;
}
// compare function for sorting
bool compare(Point p1, Point q1)
{
    Point p = Point(p1.x - mid.x,
        p1.y - mid.y);
    Point q = Point(q1.x - mid.x,
        q1.y - mid.y);
    int one = quad(p);
    int two = quad(q);
    if (one != two)
        return (one < two);
    return (p.y * q.x < q.y* p.x);
}
// Finds upper tangent of two polygons 'a' and 'b'
// represented as two vectors.
vector<Point> merger(vector<Point> a,
    vector<Point> b)
{
    // n1 -> number of points in polygon a
    // n2 -> number of points in polygon b
    int n1 = a.size(), n2 = b.size();
    int ia = 0, ib = 0;
    for (int i = 1; i < n1; i++)
        if (a[i].x > a[ia].x)
            ia = i;
    // ib -> leftmost point of b
    for (int i = 1; i < n2; i++)
        if (b[i].x < b[ib].x)
            ib = i;
    // finding the upper tangent
    int inda = ia, indb = ib;
    bool done = 0;
    while (!done)
    {
        done = 1;
        while (orientation(b[indb], a[inda], a[(inda + 1) % n1]) >= 0)
            inda = (inda + 1) % n1;
        while (orientation(a[inda], b[indb], b[(n2 + indb - 1) % n2]) <= 0)
        {
            indb = (n2 + indb - 1) % n2;
            done = 0;
        }
    }
    int uppera = inda, upperb = indb;
    inda = ia, indb = ib;
    done = 0;
    int g = 0;
    while (!done)//finding the lower tangent
    {
        done = 1;
        while (orientation(a[inda], b[indb], b[(indb + 1) % n2]) >= 0)
            indb = (indb + 1) % n2;
        while (orientation(b[indb], a[inda], a[(n1 + inda - 1) % n1]) <= 0)
        {
            inda = (n1 + inda - 1) % n1;
            done = 0;
        }
    }
    int lowera = inda, lowerb = indb;
    vector<Point> ret;
    //ret contains the convex hull after merging the two convex hulls
    //with the points sorted in anti-clockwise order
    int ind = uppera;
    ret.push_back(a[uppera]);
    while (ind != lowera)
    {
        ind = (ind + 1) % n1;
        ret.push_back(a[ind]);
    }
    ind = lowerb;
    ret.push_back(b[lowerb]);
    while (ind != upperb)
    {
        ind = (ind + 1) % n2;
        ret.push_back(b[ind]);
    }
    return ret;
}
// Brute force algorithm to find convex hull for a set
// of less than 6 points
vector<Point> bruteHull(vector<Point> a)
{
    // Take any pair of points from the set and check
    // whether it is the edge of the convex hull or not.
    // if all the remaining points are on the same side
    // of the line then the line is the edge of convex
    // hull otherwise not
    set<Point>s;
    for (int i = 0; i < a.size(); i++)
    {
        for (int j = i + 1; j < a.size(); j++)
        {
            int x1 = a[i].x, x2 = a[j].x;
            int y1 = a[i].y, y2 = a[j].y;
            int a1 = y1 - y2;
            int b1 = x2 - x1;
            int c1 = x1 * y2 - y1 * x2;
            int pos = 0, neg = 0;
            for (int k = 0; k < a.size(); k++)
            {
                if (a1 * a[k].x + b1 * a[k].y + c1 <= 0)
                    neg++;
                if (a1 * a[k].x + b1 * a[k].y + c1 >= 0)
                    pos++;
            }
            if (pos == a.size() || neg == a.size())
            {
                s.insert(a[i]);
                s.insert(a[j]);
            }
        }
    }
    vector<Point>ret;
    for (auto e : s)
        ret.push_back(e);
    // Sorting the points in the anti-clockwise order
    mid = { 0, 0 };
    int n = ret.size();
    for (int i = 0; i < n; i++)
    {
        mid.x += ret[i].x;
        mid.y += ret[i].y;
        ret[i].x *= n;
        ret[i].y *= n;
    }
    sort(ret.begin(), ret.end(), compare);
    for (int i = 0; i < n; i++)
        ret[i] = Point(ret[i].x / n, ret[i].y / n);
    return ret;
}
// Returns the convex hull for the given set of points
vector<Point> divide(vector<Point>& a)
{
   // if (a.size() <= 5)
    //    return bruteHull(a);
    
    //da se dovoljno ogranici dubina rekurzije
  //  if (a.size() <= SIZE /4)
    if (a.size() <= 10000){
        cout << "osnovni slucaj" << endl;
        return convex_hull(a); //MonotoneChain algorithm
        }

    // left contains the left half points
    // right contains the right half points
    vector<Point>left, right;
    for (int i = 0; i < a.size() / 2; i++)
        left.push_back(a[i]);
    for (int i = a.size() / 2; i < a.size(); i++)
        right.push_back(a[i]);

    vector<Point>left_hull;
    vector<Point>right_hull;

#pragma omp task shared(left_hull)
        {
            left_hull = convex_hull(left); //divide(left);
        }

#pragma omp task shared(right_hull)
        {
            right_hull = convex_hull(right); // divide(right);
        }
    
    
    // merging the convex hulls
#pragma omp taskwait
   return merger(left_hull, right_hull);
}

/*
vector<Point> divideParallel(vector<Point> a)
{
    // If the number of points is less than 6 then the
    // function uses the brute algorithm to find the
    // convex hull
    if (a.size() <= 5)
        return bruteHull(a);

    // left contains the left half points
    // right contains the right half points
    vector<Point>left, right;
    for (int i = 0; i < a.size() / 2; i++)
        left.push_back(a[i]);
    for (int i = a.size() / 2; i < a.size(); i++)
        right.push_back(a[i]);

    vector<Point>left_hull;
    vector<Point>right_hull;
    // convex hull for the left and right sets
#pragma omp task shared(left_hull)
    left_hull = divide(left);
#pragma omp task shared(right_hull)
    right_hull = divide(right);

    // merging the convex hulls
#pragma omp taskwait
    return merger(left_hull, right_hull);


}*/

vector<Point> divide_and_conquer(vector<Point> a)
{
    cout << "divide and conquer" << endl;
    std::sort(std::execution::par_unseq, a.begin(), a.end());
    cout << "divide" << endl;
    return divide(a);

}


// Driver code
int main()
{
//    vector<pair<int, int> > a;
   /* a.push_back(make_pair(0, 0));
    a.push_back(make_pair(1, -4));
    a.push_back(make_pair(-1, -5));
    a.push_back(make_pair(-5, -3));
    a.push_back(make_pair(-3, -1));
    a.push_back(make_pair(-1, -3));
    a.push_back(make_pair(-2, -2));
    a.push_back(make_pair(-1, -1));
    a.push_back(make_pair(-2, -1));
    a.push_back(make_pair(-1, 1));
    //int n = a.size();*/
    int n = 100000;

    double suma = 0;
    vector<Point> ans;
    vector<Point> ans_compare;
  /*  std::sort(a.begin(), a.end());
    ans = divide(a);*/

    cout << "enter no of threads" << endl;
    int xx;
    cin >> xx;
    omp_set_num_threads(xx);
   
    
    int x, y;
    srand(time(0));
    for (int i = 0; i < 1; i++) {
        vector<Point> a;
        vector<Point> b;
        cout << "generisanje brojeva" << endl;
        for (llu i = 0;i < SIZE;i++) {
           x = rand() % n;
           y = rand() % n;
           a.push_back(Point(x, y));
           b.push_back(Point(x, y));
        }
        cout << "kraj generisanja brojeva" << endl;
        // Find the convex hull
        double e;
        {
            Timer t("Vrijeme na n threadova: ");
            ans = divide_and_conquer(a); //divide_and_conquer(a);
            e = t.elapsed();
        }
        cout << ": " << e / 1000000 << " seconds" << endl;
        suma += e / 1000000;

        if (i == 0) {
            ans_compare = convex_hull(b);
        }
    }

  //  cout << "Prosjek: " << suma / 5 << endl;

   // vector<pair<int, int> >ans = divide_and_conquer(a);
    cout << "convex hull:\n";
    for (auto e : ans)
        cout << e.x << " "
        << e.y << endl;

    cout << "convex hull sa monotone chain:\n";
    for (auto e : ans_compare)
        cout << e.x << " "
        << e.y << endl;

    cout << "razlika:\n";
    for (auto e : ans_compare) {
        if (std::find(ans.begin(), ans.end(), e) != ans.end()) {
            /* v contains x */
        }
        else {
            cout <<"prvi skup ne sadrzi tacku "<< e.x << " "
                << e.y << endl;
        }
     
    }

    return 0;
}

