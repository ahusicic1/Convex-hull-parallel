// MyLastStrawAKAChanovAlgoritmyyy.cpp : This file
//  contains my last straw to finish this project. 
// Program execution begins and ends there
// and so does my motivation for this project.
//
#include <iostream>
#include <stack> 
#include <stdlib.h> 
#include <omp.h>
#include <time.h>
#include <chrono>
#include<vector>
#include <algorithm>
#include <execution>
#include <numeric>

using namespace std;
using namespace std::chrono;
#define ll long long int 
#define SIZE 50000
struct Point
{
	ll x, y;
};
Point p0;
ll n = SIZE;
Point points[SIZE];
Point nextToTop(stack<Point>& S)
{
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

void swap(Point& p1, Point& p2)
{
	Point temp = p1;
	p1 = p2;
	p2 = temp;
}
ll distSq(Point p1, Point p2)
{
	return (p1.x - p2.x) * (p1.x - p2.x) +
		(p1.y - p2.y) * (p1.y - p2.y);
}

// To find orientation of ordered triplet (p, q, r). 
// The function returns following values 
// 0 --> p, q and r are colinear 
// 1 --> Clockwise 
// 2 --> Counterclockwise 
ll orientation(Point p, Point q, Point r)
{
	ll val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0; // colinear 
	return (val > 0) ? 1 : 2; // clock or counterclock wise 
}

int compare(const void* vp1, const void* vp2)
{
	Point* p1 = (Point*)vp1;
	Point* p2 = (Point*)vp2;
	ll o = orientation(p0, *p1, *p2);
	if (o == 0)
		return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

	return (o == 2) ? -1 : 1;
}

// Prints convex hull of a set of n points.
stack<Point> convexHull(Point points[], int n)
{
	// Find the bottommost point
	int ymin = points[0].y, min = 0;
	for (int i = 1; i < n; i++)
	{
		int y = points[i].y;

		// Pick the bottom-most or choose the left
		// most point in case of tie
		if ((y < ymin) || (ymin == y &&
			points[i].x < points[min].x))
			ymin = points[i].y, min = i;
	}

	// Place the bottom-most point at first position
	swap(points[0], points[min]);

	// Sort n-1 points with respect to the first point.
	// A point p1 comes before p2 in sorted output if p2
	// has larger polar angle (in counterclockwise
	// direction) than p1
	p0 = points[0];
	qsort(&points[1], n - 1, sizeof(Point), compare);

	// If two or more points make same angle with p0,
	// Remove all but the one that is farthest from p0
	// Remember that, in above sorting, our criteria was
	// to keep the farthest point at the end when more than
	// one points have same angle.
	int m = 1; // Initialize size of modified array
	for (int i = 1; i < n; i++)
	{
		// Keep removing i while angle of i and i+1 is same
		// with respect to p0
		while (i < n - 1 && orientation(p0, points[i],
			points[i + 1]) == 0)
			i++;


		points[m] = points[i];
		m++;  // Update size of modified array
	}

	// If modified array of points has less than 3 points,
	// convex hull is not possible
	if (m < 3) {
		stack<Point> S;
		for (int i = 0; i < m; i++) {
			S.push(points[i]);
		}
		return S;
	}

	// Create an empty stack and push first three points
	// to it.
	stack<Point> S;
	S.push(points[0]);
	S.push(points[1]);
	S.push(points[2]);

	// Process remaining n-3 points
	for (int i = 3; i < m; i++)
	{
		// Keep removing top while the angle formed by
		// points next-to-top, top, and points[i] makes
		// a non-left turn
		while (S.size() > 1 && orientation(nextToTop(S), S.top(), points[i]) != 2)
			S.pop();
		S.push(points[i]);
	}

	return S;
}

// Driver program to test above functions 
int main()
{
	Point points[SIZE];
	int n = 100;
	for (ll i = 0;i < SIZE;i++) {
		points[i].x = rand() % n;
		points[i].y = rand() % n;

		//cout << points[i].x << ", " << points[i].y << endl;
	}
	int n1 = sizeof(points) / sizeof(points[0]);
	cout << "n1 je " << n1;
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	stack<Point> s = convexHull(points, n1);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "Elapsed time: " << time_span.count() << " seconds.";
	return 0;
	/*Point points[SIZE];
	Point points2[SIZE];
	long end_time, start_time;
	double time_overhead;
	int n = 100;
	for (ll i = 0;i < n;i++) {
		points[i].x = rand() % n;
		points[i].y = rand() % n;
	}
	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	stack<Point> s = convexHull(points, SIZE);
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
	cout << "Elapsed time: " << time_span.count() << " seconds.";

	return 0;*/
}