#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
using namespace std;


class number_type
{
public:

	number_type(void)
	{
		vertex_data.resize(vertex_length, 0);
	}

	number_type(size_t var_count)
	{
		vertex_length = var_count;
		vertex_data.resize(vertex_length, 0);
	}

	float magnitude(void)
	{
		float all_self_dot = 0;

		for (size_t i = 0; i < vertex_length; i++)
			all_self_dot += (vertex_data[i] * vertex_data[i]);

		return sqrtf(all_self_dot);
	}
	
	number_type operator+(const number_type& right) const
	{
		number_type out(right.vertex_length);

		for(size_t i = 0; i < right.vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] + right.vertex_data[i];

		return out;
	}

	number_type operator/(const float& right) const
	{
		number_type out(vertex_length);

		for (size_t i = 0; i < vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] / right;

		return out;
	}

	size_t vertex_length = 1; // default is one float
	vector<float> vertex_data;
};

number_type conj_number_type(number_type& in)
{
	number_type out;

	out.vertex_data[0] = in.vertex_data[0];

	for (size_t i = 1; i < in.vertex_length; i++)
		out.vertex_data[i] = -in.vertex_data[i];

	return out;
}

number_type pow_number_type(number_type &in, float exponent)
{
	const float beta = exponent;

	const float fabs_beta = fabsf(beta);

	float all_self_dot = 0;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	float imag_self_dot = 0;

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	number_type out;

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float all_len = sqrtf(all_self_dot);
	const float imag_len = sqrtf(imag_self_dot);
	const float		 self_dot_beta = powf(all_self_dot, fabs_beta / 2.0f);

	out.vertex_data[0] = self_dot_beta * std::cos(fabs_beta * std::acos(in.vertex_data[0] / all_len));

	for (size_t i = 1; i < out.vertex_length; i++)
		out.vertex_data[i] = in.vertex_data[i] * self_dot_beta * sin(fabs_beta * acos(in.vertex_data[0] / all_len)) / imag_len;

	if (beta < 0)
		out = conj_number_type(out) / powf(out.magnitude(), 2.0f);

	return out;
}

inline float iterate(
	number_type Z,
	const number_type C,
	const short unsigned int max_iterations,
	const float threshold)
{
	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = pow_number_type(Z, 2.0) + C;

		if (Z.magnitude() >= threshold)
			break;
	}

	return Z.magnitude();
}


int main(void)
{
	const size_t res = 10;
	const float grid_max = 1.5;
	const float grid_min = -grid_max;
	const unsigned short int max_iterations = 8;
	const float threshold = 4.0f;
	const float step_size = (grid_max - grid_min) / (res - 1);

	number_type C(5);

	srand(1);

	for (size_t i = 0; i < C.vertex_length; i++)
		C.vertex_data[i] = rand() / static_cast<float>(RAND_MAX) * 0.5f;

	number_type Z(5);

	for (size_t i = 0; i < Z.vertex_length; i++)
		Z.vertex_data[i] = grid_min;

	size_t total_count = 0;
	size_t in_set = 0;

	for (size_t i0 = 0; i0 < res; i0++, Z.vertex_data[0] += step_size)
	{
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
							in_set++;
						// use z here
					}
				}
			}
		}
	}

	cout << in_set << " of " << total_count << endl;






	return 0;
}