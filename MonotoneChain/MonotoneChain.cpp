// C++ implementation of the approach
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
#include <execution>

using namespace std;
using namespace std::chrono;
#define llu long long int
#define llu long long int
using namespace std;

#define SIZE 50000000

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

struct Point {

	llu x, y;

	bool operator<(Point p)
	{
		return x < p.x || (x == p.x && y < p.y);
	}
};

// Cross product of two vectors OA and OB
// returns positive for counter clockwise
// turn and negative for clockwise turn
llu cross_product(Point O, Point A, Point B)
{
	return (A.x - O.x) * (B.y - O.y)
		- (A.y - O.y) * (B.x - O.x);
}

// Returns a list of points on the convex hull
// in counter-clockwise order
vector<Point> convex_hull(vector<Point>& A)
{

	cout << "Convex hull monotone chain" << endl;
	int n = A.size(), k = 0;

	if (n <= 3)
		return A;

	vector<Point> ans(2 * n);

	// Sort points lexicographically
	//std::sort(A.begin(), A.end());
	std::sort(std::execution::par_unseq, A.begin(), A.end());


	// Build lower hull
	for (int i = 0; i < n; ++i) {

		// If the point at K-1 position is not a part
		// of hull as vector from ans[k-2] to ans[k-1]
		// and ans[k-2] to A[i] has a clockwise turn
		{
			while (k >= 2 && cross_product(ans[k - 2],
				ans[k - 1], A[i]) <= 0)
				k--;
			ans[k++] = A[i];
		}
	}

	// Build upper hull
	for (size_t i = n - 1, t = k + 1; i > 0; --i) {

		// If the point at K-1 position is not a part
		// of hull as vector from ans[k-2] to ans[k-1]
		// and ans[k-2] to A[i] has a clockwise turn
		while (k >= t && cross_product(ans[k - 2],
			ans[k - 1], A[i - 1]) <= 0)
			k--;
		ans[k++] = A[i - 1];
	}

	// Resize the array to desired size
	ans.resize(k - 1);

	return ans;
}


// Driver code
int main()
{
	//vector<Point> points;
	/*
	// Add points
	points.push_back({ 0, 3 });
	points.push_back({ 2, 2 });
	points.push_back({ 1, 1 });
	points.push_back({ 2, 1 });
	points.push_back({ 3, 0 });
	points.push_back({ 0, 0 });
	points.push_back({ 3, 3 });*/

	vector<Point> points(SIZE);

	//vece n znaci duze vrijeme generisanja tacaka (koje nam svakako nije bitno)
	//ali je bitno da je n dovoljno veliko u odnosu na velicinu points kako
	//bi convex hull bio vise random 
	//za n << SIZE gotovo je sigurno da ce se generisati taèke koje ogranièavaju kvadrat
	//dopustivog prostora tako da æe u tom sluèaju samo te 4 taèke èiniti convex hull
	int n = 100000; 

	double suma = 0;
	vector<Point> ans;

	cout << "enter no of threads" << endl;
	int xx;
	cin >> xx;
	omp_set_num_threads(xx);

	srand(time(0));
	for (int i = 0; i < 5; i++) {
		for (llu i = 0;i < SIZE;i++) {
			points[i].x = rand() % n;
			points[i].y = rand() % n;
		}
		
		// Find the convex hull
		double e;
		{
			Timer t("Vrijeme na n threadova: ");
			ans = convex_hull(points);
			 e = t.elapsed();
		}
		cout << ": " << e / 1000000 << " seconds" << endl;
		suma += e / 1000000;
	}

	cout << "Prosjek: " << suma / 5;

	// Print the convex hull
	for (int i = 0; i < ans.size(); i++)
		cout << "(" << ans[i].x << ", "
		<< ans[i].y << ")" << endl;


	return 0;
}
