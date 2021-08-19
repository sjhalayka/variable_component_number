#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
using namespace std;


class sedenion
{
public:

	sedenion(void)
	{
		vertex_data.resize(vertex_length, 0);
	}

	size_t vertex_length = 16;
	vector<float> vertex_data;
};


sedenion pow_sedenion(sedenion &in, float exponent)
{
	const float beta = exponent;

	const float fabs_beta = fabsf(beta);

	float all_self_dot = 0;

	for (size_t i = 0; i < in.vertex_length; i++)
		all_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	float imag_self_dot = 0;

	for (size_t i = 1; i < in.vertex_length; i++)
		imag_self_dot += (in.vertex_data[i] * in.vertex_data[i]);

	sedenion out;

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

	return out;

	//if (beta < 0)
	//	inverse(&out, 0, &out);
}


complex<float> pow_complex(const complex<float>& in, const float beta)
{
	float fabs_beta = fabsf(beta);

	float self_dot = in.real() * in.real() + in.imag() * in.imag();

	if (self_dot == 0)
		return complex<float>(0, 0);

	float len = std::sqrtf(self_dot);
	float self_dot_beta = std::powf(self_dot, fabs_beta / 2.0f);

	complex<float> out = complex<float>(
		self_dot_beta * std::cos(fabs_beta * std::acos(in.real() / len)),
		in.imag() * self_dot_beta * std::sin(fabs_beta * std::acos(in.real() / len)) / sqrtf(in.imag() * in.imag()));

	//if (beta < 0)
	//	out = conj(out) / powf(abs(out), 2.0f);

	return out;

}




int main(void)
{
	sedenion s;

	s.vertex_data[0] = 0.2f;
	s.vertex_data[1] = 0.5f;
	s = pow_sedenion(s, 2.0f);

	cout << s.vertex_data[0] << " " << s.vertex_data[1] << endl;


	complex<float> c(0.2f, 0.5f);
	c = pow_complex(c, 2.0);

	cout << c.real() << " " << c.imag() << endl;

	return 0;
}