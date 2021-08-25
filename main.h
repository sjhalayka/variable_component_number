#pragma once


class vertex_4
{
public:
	float x, y, z, w;
};

class vertex_3
{
public:

	inline vertex_3(void) : x(0.0f), y(0.0f), z(0.0f) { /*default constructor*/ }
	inline vertex_3(const float src_x, const float src_y, const float src_z, const size_t src_index) : x(src_x), y(src_y), z(src_z) { /* custom constructor */ }

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
};



class quintonion
{
public:

	quintonion(void)
	{
		vertex_data.resize(vertex_length, 0);
	}

	float magnitude(void)
	{
		float all_self_dot = 0;

		for (size_t i = 0; i < vertex_length; i++)
			all_self_dot += (vertex_data[i] * vertex_data[i]);

		return sqrtf(all_self_dot);
	}

	quintonion operator+(const quintonion& right) const
	{
		quintonion out;

		for (size_t i = 0; i < right.vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] + right.vertex_data[i];

		return out;
	}

	quintonion operator/(const float& right) const
	{
		quintonion out;

		for (size_t i = 0; i < vertex_length; i++)
			out.vertex_data[i] = vertex_data[i] / right;

		return out;
	}

	size_t vertex_length = 5;
	vector<float> vertex_data;
};

quintonion conj_number_type(quintonion& in)
{
	quintonion out;

	out.vertex_data[0] = in.vertex_data[0];

	for (size_t i = 1; i < in.vertex_length; i++)
		out.vertex_data[i] = -in.vertex_data[i];

	return out;
}

quintonion pow_number_type(quintonion& in, float exponent)
{
	const float beta = exponent;
	const float fabs_beta = fabsf(beta);
	float all_self_dot = 0;
	float imag_self_dot = 0;
	quintonion out;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	if (all_self_dot == 0)
	{
		for (size_t i = 0; i < out.vertex_length; i++)
			out.vertex_data[i] = 0;

		return out;
	}

	const float all_len = sqrtf(all_self_dot);
	const float imag_len = sqrtf(imag_self_dot);
	const float self_dot_beta = powf(all_self_dot, fabs_beta / 2.0f);

	out.vertex_data[0] = self_dot_beta * std::cos(fabs_beta * std::acos(in.vertex_data[0] / all_len));

	for (size_t i = 1; i < out.vertex_length; i++)
		out.vertex_data[i] = in.vertex_data[i] * self_dot_beta * sin(fabs_beta * acos(in.vertex_data[0] / all_len)) / imag_len;

	if (beta < 0)
		out = conj_number_type(out) / powf(out.magnitude(), 2.0f);

	return out;
}

inline float iterate(
	quintonion Z,
	quintonion C,
	const short unsigned int max_iterations,
	const float threshold)
{
	//C = Z;

	//Z.vertex_data[0] = 0.0f;
	//Z.vertex_data[1] = 0.0f;
	//Z.vertex_data[2] = 0.0f;
	//Z.vertex_data[3] = 0.0f;
	//Z.vertex_data[4] = 0.0f;


	for (short unsigned int i = 0; i < max_iterations; i++)
	{
		Z = pow_number_type(Z, 2.0) + C;

		if (Z.magnitude() >= threshold)
			break;
	}

	return Z.magnitude();
}