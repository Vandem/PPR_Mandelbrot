// disable visual studio warning
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <chrono>
#include "bitmap_image.hpp"

using namespace std;

std::tuple<int, int, int> HSVtoRGB(float H, float S, float V) {
	if (H > 360 || H < 0 || S>100 || S < 0 || V>100 || V < 0) {
		cout << "The givem HSV values are not in valid range" << endl;
		return {0, 0, 0};
	}
	float s = S / 100;
	float v = V / 100;
	float C = s * v;
	float X = C * (1 - abs(fmod(H / 60.0, 2) - 1));
	float m = v - C;
	float r, g, b;
	if (H >= 0 && H < 60) {
		r = C, g = X, b = 0;
	}
	else if (H >= 60 && H < 120) {
		r = X, g = C, b = 0;
	}
	else if (H >= 120 && H < 180) {
		r = 0, g = C, b = X;
	}
	else if (H >= 180 && H < 240) {
		r = 0, g = X, b = C;
	}
	else if (H >= 240 && H < 300) {
		r = X, g = 0, b = C;
	}
	else {
		r = C, g = 0, b = X;
	}

	int R = (r + m) * 255;
	int G = (g + m) * 255;
	int B = (b + m) * 255;

	return { R, G, B };
}

std::tuple<double, double> normalizeToViewRect(int px, int py, double dx, double dy, int minX, int minY) {
	double x = minX + px * dx;
	double y = minY + py * dy;
	return { x, y };
}

void calculatePixel2(int xImage, int yImage, double cx, double cy, int max_iterations, bitmap_image& image) {
	double zx = cx;
	double zy = cy;

	for (unsigned int i = 0; i < max_iterations; ++i)
	{
		double x = zx * zx - zy * zy + cx;
		double y = 2 * zx * zy + cy;

		if (((x * x) + (y * y)) > 4)
		{
			const double z = sqrt(x * x + y * y);

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

void mandelbrot2() {
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
	cout << "Computing...";


	bitmap_image image(w, h);
	image.clear();

	auto start = chrono::high_resolution_clock::now();

	double dx = (maxX - minX) / ((double)image.width() - 1);
	double dy = (maxY - minY) / ((double)image.height() - 1);

	int z, yImage, xImage, i;
	double cx, cy, zx, x, zy, y;

# pragma omp parallel \
  shared ( max_iterations, maxX, minX, maxY, minY, z, image ) \
  private ( yImage, xImage, i, cx, cy, zx, x, zy, y )
	{
# pragma omp for
		for (yImage = 0; yImage < image.height(); ++yImage)
		{
			for (xImage = 0; xImage < image.width(); ++xImage)
			{
				auto [cx, cy] = normalizeToViewRect(xImage, yImage, dx, dy, minX, minY);
				//calculatePixel2(xImage, yImage, cx, cy, max_iterations, image);

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

	image.save_image("mandelbrot_set.bmp");
}

int calculatePixel1(double real, double imag) {
	int limit = 100;
	double zReal = real;
	double zImag = imag;

	for (int i = 0; i < limit; ++i) {
		double r2 = zReal * zReal;
		double i2 = zImag * zImag;

		if (r2 + i2 > 4.0) return i;

		zImag = 2.0 * zReal * zImag + imag;
		zReal = r2 - i2 + real;
	}
	return limit;
}

void mandelbrot1() {
	int width = 1024;
	int height = 1024;

	double minX = -2.0;
	double minY = -1.0;
	double maxX = 1.0;
	double maxY = 1.0;

	double dx = (maxX - minX) / (width - 1);
	double dy = (maxY - minY) / (height - 1);

	FILE* f = fopen("out.ppm", "wb");
	fprintf(f, "P6\n%i %i 255\n", width, height);

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {

			double x = minX + j * dx;
			double y = maxY - i * dy;

			int value = calculatePixel1(x, y);

			int s = 100, v = 100;
			if (value == 100) { s = 0, v = 0; };
			int scaledValue = value * (255 / 100);
			auto [r, g, b] = HSVtoRGB(scaledValue, s, v);

			fputc(r, f);
			fputc(g, f);
			fputc(b, f);

		}
	}

	fclose(f);
}


int main() {

	//mandelbrot1();
	mandelbrot2();
	return 0;
}
