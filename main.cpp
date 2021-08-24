#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
using namespace std;

#include "main.h"







int main(void)
{
	ofstream out_bin("points.bin", ios_base::binary);

	const size_t res = 25;
	const float grid_max = 1.5;
	const float grid_min = -grid_max;
	const unsigned short int max_iterations = 8;
	const float threshold = 4.0f;
	const float step_size = (grid_max - grid_min) / (res - 1);

	quintonion C;
	C.vertex_data[0] = 0.2f;
	C.vertex_data[1] = 0.5f;
	C.vertex_data[2] = 0.3f;
	C.vertex_data[3] = 0.2f;
	C.vertex_data[4] = 0.1f;

	quintonion Z;

	for (size_t i = 0; i < Z.vertex_length; i++)
		Z.vertex_data[i] = grid_min;

	size_t total_count = 0;
	size_t in_set = 0;

	for (size_t i0 = 0; i0 < res; i0++, Z.vertex_data[0] += step_size)
	{
		cout << i0 + 1 << " of " << res << endl;

		Z.vertex_data[1] = grid_min;

		for (size_t i1 = 0; i1 < res; i1++, Z.vertex_data[1] += step_size)
		{
			Z.vertex_data[2] = grid_min;

			for (size_t i2 = 0; i2 < res; i2++, Z.vertex_data[2] += step_size)
			{
				Z.vertex_data[3] = grid_min;

				for (size_t i3 = 0; i3 < res; i3++, Z.vertex_data[3] += step_size)
				{
					Z.vertex_data[4] = grid_min;

					for (size_t i4 = 0; i4 < res; i4++, Z.vertex_data[4] += step_size)
					{
						total_count++;

						float magnitude = iterate(Z, C, max_iterations, threshold);

						if (magnitude < threshold)
						{
							in_set++;

							vertex_4 v4;

							if (0)//Z.vertex_data[4] != 0)
							{
								v4.x = Z.vertex_data[0] / Z.vertex_data[4];
								v4.y = Z.vertex_data[1] / Z.vertex_data[4];
								v4.z = Z.vertex_data[2] / Z.vertex_data[4];
								v4.w = Z.vertex_data[3] / Z.vertex_data[4];
							}
							else
							{
								v4.x = Z.vertex_data[0];
								v4.y = Z.vertex_data[1];
								v4.z = Z.vertex_data[2];
								v4.w = Z.vertex_data[3];
							}

							vertex_3 v;

							if (0)//v4.w != 0)
							{
								v.x = v4.x / v4.w;
								v.y = v4.y / v4.w;
								v.z = v4.z / v4.w;
							}
							else
							{
								v.x = v4.x;
								v.y = v4.y;
								v.z = v4.z;
							}

							out_bin.write(reinterpret_cast<const char*>(&v.x), sizeof(float));
							out_bin.write(reinterpret_cast<const char*>(&v.y), sizeof(float));
							out_bin.write(reinterpret_cast<const char*>(&v.z), sizeof(float));
						}
					}
				}
			}
		}
	}

	out_bin.close();

	cout << in_set << " of " << total_count << endl;

	return 0;
}