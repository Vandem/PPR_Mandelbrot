// disable visual studio warning
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <chrono>
#include <omp.h>
#include <vector>
#include <numeric>
#include <execution>
#include "bitmap_image.hpp"

#define NUM_THREADS 8

using namespace std;

std::tuple<double, double> normalizeToViewRect(int px, int py, double dx, double dy, int minX, int minY) {
	double x = minX + px * dx;
	double y = minY + py * dy;
	return { x, y };
}

long mandelbrot(bool saveImage = true) {
	unsigned int w = 1024, h = 1024;
	int minX = -2;
	int minY = -1;
	int maxX = 1;
	int maxY = 1;
	unsigned int max_iterations = 1000;

	cout << "Mandelbrot set. Start by providing the necessary parameters.\n";
	cout << "Image width: ";
	cin >> w;
	cout << "Image height: ";
	cin >> h;
	cout << "minX: ";
	cin >> minX;
	cout << "minY: ";
	cin >> minY;
	cout << "maxX: ";
	cin >> maxX;
	cout << "maxY: ";
	cin >> maxY;
	cout << "Max iterations: ";
	cin >> max_iterations;
	cout << "Computing...\n";


	bitmap_image image(w, h);
	image.clear();

	int z, yImage, xImage, i;
	double cx, cy, zx, x, zy, y, dx, dy;

	auto start = chrono::high_resolution_clock::now();

	dx = (maxX - minX) / ((double)w - 1);
	dy = (maxY - minY) / ((double)h - 1);



	omp_set_num_threads(NUM_THREADS);
# pragma omp parallel \
  shared ( max_iterations, maxX, minX, maxY, minY, w, h) \
  private ( yImage, xImage, i, cx, cy, zx, x, zy, y, z )
	{
# pragma omp for
		for (yImage = 0; yImage < h; ++yImage)
		{
			for (xImage = 0; xImage < w; ++xImage)
			{
				auto [cx, cy] = normalizeToViewRect(xImage, yImage, dx, dy, minX, minY);

				zx = cx;
				zy = cy;

				for (int i = 0; i < max_iterations; ++i)
				{
					x = zx * zx - zy * zy + cx;
					y = 2 * zx * zy + cy;

					if (((x * x) + (y * y)) > 4)
					{
						z = sqrt(x * x + y * y);

						//https://en.wikipedia.org/wiki/Mandelbrot_set#Continuous_.28smooth.29_coloring
						const unsigned int index = static_cast<unsigned int>
							(1000.0 * log2(1.75 + i - log2(log2(z))) / log2(max_iterations));

						image.set_pixel(xImage, yImage, jet_colormap[index]);

						break;
					}
					zx = x;
					zy = y;
				}
			}
		}
	}

	auto stop = chrono::high_resolution_clock::now();
	auto duration = chrono::duration_cast<chrono::milliseconds>(stop - start);
	cout << "Time taken by function: "
		<< duration.count() << " milliseconds" << endl;

	if (saveImage) {
		image.save_image("mandelbrot_set.bmp");
	}

	return duration.count();
}

long getAverage(std::vector<long> const& v){
    if(v.empty()){
        return 0;
    }

    auto const count = static_cast<long>(v.size());
    return reduce(v.begin(), v.end()) / count;
}


int main() {

	mandelbrot();

	// code used for benchmarking
	//vector<long> durations;
	//for (int i = 0; i < 100; i++)
	//{
	//	durations.push_back(mandelbrot(false));
	//}
	//long average = getAverage(durations);
	//printf("%d", average);
	
	return 0;
}
