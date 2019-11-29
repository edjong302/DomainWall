#ifndef INITIAL_CONDITIONS
#define INITIAL_CONDITIONS

class initial_conditions{
	private:
		float offset1 = 0.25;
		float offset2 = 0.75;
		float sigma = 0.05;
	public:
		float get_offset1();
		float get_offset2();
		float get_sigma();
		float phi_init(float i_offset1, float i_offset2, float i_sigma, float i_x);
		float psi_init(float i_offset1, float i_offset2, float i_sigma, float i_x);
};

#endif