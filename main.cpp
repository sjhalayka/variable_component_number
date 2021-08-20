#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
using namespace std;

class vertex_3
{
public:

	inline vertex_3(void) : x(0.0f), y(0.0f), z(0.0f), index(0) { /*default constructor*/ }
	inline vertex_3(const float src_x, const float src_y, const float src_z, const size_t src_index) : x(src_x), y(src_y), z(src_z), index(src_index) { /* custom constructor */ }

	inline bool operator==(const vertex_3& right) const
	{
		if (right.x == x && right.y == y && right.z == z)
			return true;
		else
			return false;
	}

	inline bool operator<(const vertex_3& right) const
	{
		if (right.x > x)
			return true;
		else if (right.x < x)
			return false;

		if (right.y > y)
			return true;
		else if (right.y < y)
			return false;

		if (right.z > z)
			return true;
		else if (right.z < z)
			return false;

		return false;
	}

	inline const vertex_3& operator-(const vertex_3& right) const
	{
		static vertex_3 temp;

		temp.x = this->x - right.x;
		temp.y = this->y - right.y;
		temp.z = this->z - right.z;

		return temp;
	}

	inline const vertex_3& operator+(const vertex_3& right) const
	{
		static vertex_3 temp;

		temp.x = this->x + right.x;
		temp.y = this->y + right.y;
		temp.z = this->z + right.z;

		return temp;
	}

	inline const vertex_3& operator*(const float& right) const
	{
		static vertex_3 temp;

		temp.x = this->x * right;
		temp.y = this->y * right;
		temp.z = this->z * right;

		return temp;
	}

	inline const vertex_3& cross(const vertex_3& right) const
	{
		static vertex_3 temp;

		temp.x = y * right.z - z * right.y;
		temp.y = z * right.x - x * right.z;
		temp.z = x * right.y - y * right.x;

		return temp;
	}

	inline float dot(const vertex_3& right) const
	{
		return x * right.x + y * right.y + z * right.z;
	}

	inline const float self_dot(void)
	{
		return x * x + y * y + z * z;
	}

	inline const float length(void)
	{
		return std::sqrtf(self_dot());
	}

	inline const void normalize(void)
	{
		float len = length();

		if (0.0f != len)
		{
			x /= len;
			y /= len;
			z /= len;
		}
	}

	float x, y, z;
	size_t index;
};



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

		for (size_t i = 0; i < right.vertex_length; i++)
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

number_type pow_number_type(number_type& in, float exponent)
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
	//srand(123);// time(0));

	const size_t res = 25;
	const float grid_max = 1.5;
	const float grid_min = -grid_max;
	const unsigned short int max_iterations = 8;
	const float threshold = 4.0f;
	const float step_size = (grid_max - grid_min) / (res - 1);

	number_type C(5);
	C.vertex_data[0] = 0.2f;
	C.vertex_data[1] = 0.5f;
	C.vertex_data[2] = 0.3f;
	C.vertex_data[3] = 0.2f;
	C.vertex_data[4] = 0.1f;

	//for (size_t i = 0; i < C.vertex_length; i++)
	//	C.vertex_data[i] = rand() / static_cast<float>(RAND_MAX) * 0.5f;

	number_type Z(5);

	for (size_t i = 0; i < Z.vertex_length; i++)
		Z.vertex_data[i] = grid_min;

	size_t total_count = 0;
	size_t in_set = 0;

	set<vertex_3> all_points;

	for (size_t i0 = 0; i0 < res; i0++, Z.vertex_data[0] += step_size)
	{
		cout << i0 << " of " << res << endl;

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

							// Simply drop the 4th and 5th dimension
							vertex_3 v;
							v.x = Z.vertex_data[0];
							v.y = Z.vertex_data[1];
							v.z = Z.vertex_data[2];

							all_points.insert(v);
						}
					}
				}
			}
		}
	}

	if (all_points.size() > 0)
	{
		ofstream outfile("points.txt");

		for (set<vertex_3>::const_iterator ci = all_points.begin(); ci != all_points.end(); ci++)
			outfile << ci->x << " " << ci->y << " " << ci->z << endl;

		outfile.close();
	}
	else
	{
		cout << "Empty set" << endl;
	}
//	cout << in_set << " of " << total_count << endl;






	return 0;
}